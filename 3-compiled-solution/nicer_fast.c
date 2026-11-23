/* nicer_fast.c: one likelihood working-out of the public quick example, compiled.
 *
 * A line-by-line port of our own Python programs (2-python-solution/code/fast_kernel.py: spot_counts with the one fit of light's
 * path, and public_likelihood.py: the public likelihood's definition). Nothing of the public code is in here.
 * The tables it reads (the fit of light's path made from our own integrals, Gauss and Chebyshev tools, the energy matrix, the
 * example's data) are written once by tools/export_for_c.py.
 *
 * Light's path comes from the clock rate, the rulers, light at c in each frame at rest and the constant b along a ray
 * (see code/bend.py for the derivation). The star's flattened shape is the public example's input.
 *
 * usage: nicer_fast <tables> <rings per panel> <points> <energies> <repeats> mass radius distance cos_i phase_shift colat rho log10T [dump file]
 *        nicer_fast <tables> batch <rings per panel> <points> <energies> <threads> <stars file> [dump file]
 * build: gcc -O3 -march=native -ffast-math -fopenmp -o nicer_fast nicer_fast.c -lm      (the Makefile builds every form)
 *        (the program holds no infinity and no not-a-number, so the compiler's fast arithmetic is safe for it)
 * A negative number of repeats times that many workings-out with a different star at each one (as a sampler asks).
 * The batch mode works out many stars at once, one star per thread at any moment; nothing inside one star runs in parallel.
 * Built with -DNICER_LIBRARY -fPIC -shared it is a library with the entry points nicer_load, nicer_choose and nicer_batch.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)
#define HUGE_NEG (-1.0e300)      /* stands for "minus infinity" (the public code's failure value is near -1e90) */
#define HUGE_POS (1.0e300)

/* ---------------------------------------------------------------- tables */
#define MAXARR 128
typedef struct { char name[48]; int nd; int dims[3]; double *v; } Arr;
static Arr arrs[MAXARR];
static int narr = 0;

static void load_tables(const char *path)
{
    FILE *fh = fopen(path, "rb");
    int count, len, i, k;
    if (!fh) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    if (fread(&count, 4, 1, fh) != 1 || count > MAXARR) { fprintf(stderr, "bad table file\n"); exit(2); }
    for (i = 0; i < count; i++) {
        Arr *a = &arrs[narr++];
        size_t total = 1;
        if (fread(&len, 4, 1, fh) != 1 || len >= 48) { fprintf(stderr, "bad name\n"); exit(2); }
        if (fread(a->name, 1, (size_t)len, fh) != (size_t)len) exit(2);
        a->name[len] = 0;
        if (fread(&a->nd, 4, 1, fh) != 1 || a->nd > 3) exit(2);
        for (k = 0; k < a->nd; k++) { if (fread(&a->dims[k], 4, 1, fh) != 1) exit(2); total *= (size_t)a->dims[k]; }
        a->v = (double *)malloc(total * sizeof(double));
        if (fread(a->v, sizeof(double), total, fh) != total) { fprintf(stderr, "short read in %s\n", a->name); exit(2); }
    }
    fclose(fh);
}

static const Arr *get(const char *name)
{
    int i;
    for (i = 0; i < narr; i++) if (strcmp(arrs[i].name, name) == 0) return &arrs[i];
    fprintf(stderr, "table %s is missing\n", name);
    exit(2);
}

static const Arr *getn(const char *stem, int n)
{
    char buf[48];
    snprintf(buf, sizeof buf, "%s_%d", stem, n);
    return get(buf);
}

#define MAXR 512          /* rings */
#define MAXN 32           /* points per ring and half-turn */
#define MAXE 48           /* energies */
#define MAXB 64           /* phase bins */
#define MAXZ 96           /* terms of the path series */
#define MAXD 257          /* dense grid */
#define NCHMAX 512

static int NU, NZ, ND, NVIS, NPP, NN, NE, NB, NCH;
static const double *FIT_CT, *FIT_CA, *VUG, *VT0, *VT0X, *VT0XX, *GLG, *GLW, *CCOS, *CSIN, *CG, *CQ0, *CTOD, *EN_E, *EN_W;
static const double *COUNTS, *SUPPORT, *FACT;
static double FIT_U0, FIT_U1, ZM, FIT_XMIN, FIT_TILT;
static double C_CM, XGM, XKPC, XKB, XBB, T_EXP, FREQ;

static void choose(int npp, int n, int nE)
{
    const Arr *a = get("fit_ct"), *m = get("fit_meta"), *c = get("consts"), *d = get("data_counts");
    NU = a->dims[0]; NZ = a->dims[1];
    FIT_CT = a->v; FIT_CA = get("fit_ca")->v;
    FIT_U0 = m->v[0]; FIT_U1 = m->v[1]; ZM = m->v[2]; FIT_XMIN = m->v[3]; FIT_TILT = m->v[4];
    a = get("vis_ug"); NVIS = a->dims[0]; VUG = a->v;
    VT0 = get("vis_t0")->v; VT0X = get("vis_t0x")->v; VT0XX = get("vis_t0xx")->v;
    NPP = npp; GLG = getn("gl_g", npp)->v; GLW = getn("gl_w", npp)->v;
    NN = n; CCOS = getn("ch_cos", n)->v; CSIN = getn("ch_sin", n)->v; CG = getn("ch_G", n)->v; CQ0 = getn("ch_q0", n)->v;
    CTOD = getn("ch_tod", n)->v; ND = (int)get("ch_meta")->v[0];
    NE = nE; EN_E = getn("en_E", nE)->v; EN_W = getn("en_W", nE)->v;
    NCH = d->dims[0]; NB = d->dims[1];
    COUNTS = d->v; SUPPORT = get("data_support")->v; FACT = get("data_fact")->v;
    C_CM = c->v[0]; XGM = c->v[1]; XKPC = c->v[2]; XKB = c->v[3]; XBB = c->v[4]; T_EXP = c->v[5]; FREQ = c->v[6];
    if (NZ > MAXZ || ND > MAXD || NN > MAXN || NE > MAXE || NB > MAXB || NCH > NCHMAX || NU > 64) { fprintf(stderr, "a table is larger than the program allows\n"); exit(2); }
}

/* ---------------------------------------------------------------- small tools */
static inline double clip(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }

static inline double cheb_sum(const double *c, int n, double x)
{
    double b1 = 0.0, b2 = 0.0, x2 = 2.0 * x, t;
    int k;
    for (k = n - 1; k >= 1; k--) { t = x2 * b1 - b2 + c[k]; b2 = b1; b1 = t; }
    return x * b1 - b2 + c[0];
}

/* the swept angle and the delay (units R/c) from a ring's own series; x = cos(alpha) */
static inline double swept(const double *cth, double x)
{
    double z = sqrt(x < 1.0 ? 1.0 - x : 0.0);
    return z * cheb_sum(cth, NZ, clip(2.0 * z / ZM - 1.0, -1.0, 1.0));
}

static inline void swept_delay(const double *cth, const double *cta, double x, double *th, double *ta)
{
    double z = sqrt(x < 1.0 ? 1.0 - x : 0.0), zeta = clip(2.0 * z / ZM - 1.0, -1.0, 1.0), x2 = 2.0 * zeta;
    double p1 = 0.0, p2 = 0.0, q1 = 0.0, q2 = 0.0, t;
    int k;
    for (k = NZ - 1; k >= 1; k--) {
        t = x2 * p1 - p2 + cth[k]; p2 = p1; p1 = t;
        t = x2 * q1 - q2 + cta[k]; q2 = q1; q1 = t;
    }
    *th = z * (zeta * p1 - p2 + cth[0]);
    *ta = z * z * (zeta * q1 - q2 + cta[0]);
}

static inline double mono(double th) { double c = cos(th); return th <= M_PI ? 1.0 - c : 3.0 + c; }

static void ring_series(double u, double *cth, double *cta)
{
    double xi = clip((2.0 * u - FIT_U0 - FIT_U1) / (FIT_U1 - FIT_U0), -1.0, 1.0), Tu[64];
    int k, z;
    Tu[0] = 1.0; Tu[1] = xi;
    for (k = 2; k < NU; k++) Tu[k] = 2.0 * xi * Tu[k - 1] - Tu[k - 2];
    for (z = 0; z < NZ; z++) { cth[z] = 0.0; cta[z] = 0.0; }
    for (k = 0; k < NU; k++) {
        const double *rt = FIT_CT + (size_t)k * NZ, *ra = FIT_CA + (size_t)k * NZ;
        double tk = Tu[k];
        for (z = 0; z < NZ; z++) { cth[z] += tk * rt[z]; cta[z] += tk * ra[z]; }
    }
}

/* ---------------------------------------------------------------- latitudes */
static inline double tab_u(double u, const double *tab)      /* linear interpolation on the even grid of u: 0.12, 0.125, ..., 0.62 */
{
    double t = (u - 0.12) * 200.0;
    int i = (int)t;
    if (i < 0) i = 0;
    if (i > 99) i = 99;
    return tab[i] + (t - i) * (tab[i + 1] - tab[i]);
}

typedef struct { double L, Req, o2, incl, sx; int kind, exact; } BrkCtx;

static inline double lat_of(const BrkCtx *c, double T)
{
    return c->kind == 0 ? T - c->incl : c->kind == 1 ? TWO_PI - c->incl - T : c->kind == 2 ? c->incl + T : c->incl - T;
}

/* a contact latitude solves lat_of(Theta(u(th), x(th))) = th, with x the edge-on direction of the tilted surface there */
static double brk_resid(double th, const BrkCtx *c)
{
    double cs = cos(th), sn = sin(th), R = c->Req * (1.0 + c->o2 * cs * cs), u = c->L / R;
    double h = -2.0 * c->Req * c->o2 * cs * sn / (R * sqrt(1.0 - u)), x = c->sx * h / sqrt(1.0 + h * h), Th;
    if (c->exact) {
        double cth[MAXZ], cta[MAXZ];
        ring_series(u, cth, cta);
        Th = swept(cth, x);
    } else Th = tab_u(u, VT0) + tab_u(u, VT0X) * x + 0.5 * tab_u(u, VT0XX) * x * x;
    return lat_of(c, Th) - th;
}

/* Latitudes at which a ring's visibility changes kind (see fast_kernel.visibility_breaks): secant steps with the swept angle
 * near x = 0 from the tables, then, where the tilt is large, further secant steps on the exact function. */
static int visibility_breaks(double L, double Req, double o2, double incl, int exact, double *lat, int *left, int *right)
{
    BrkCtx c;
    int kind, it, nbk = 0;
    c.L = L; c.Req = Req; c.o2 = o2; c.incl = incl;
    for (kind = 0; kind < 4; kind++) {
        double lo_ok, hi_ok, a, b, fa, fb, nx, den;
        int fl, fr;
        c.kind = kind; c.sx = (kind == 0 || kind == 2) ? -1.0 : 1.0; c.exact = 0;
        if (kind == 0) { lo_ok = 0.0; hi_ok = M_PI - incl; fl = 0; fr = 1; }
        else if (kind == 1) { lo_ok = M_PI - incl; hi_ok = M_PI; fl = 1; fr = 0; }
        else if (kind == 2) { lo_ok = incl; hi_ok = M_PI; fl = 1; fr = 0; }
        else { lo_ok = 0.0; hi_ok = incl; fl = 0; fr = 1; }
        lo_ok = fmax(lo_ok, 0.0) + 1e-4;
        hi_ok = fmin(hi_ok, M_PI) - 1e-4;
        if (hi_ok <= lo_ok) continue;
        a = clip(lat_of(&c, tab_u(L / Req, VT0)), lo_ok, hi_ok);
        b = a + 0.02 < hi_ok ? a + 0.02 : a - 0.02;
        fa = brk_resid(a, &c); fb = brk_resid(b, &c);
        for (it = 0; it < 12; it++) {
            if (fabs(fb) < 1e-11 || fb == fa) break;
            nx = clip(b - fb * (b - a) / (fb - fa), lo_ok, hi_ok);
            a = b; fa = fb;
            b = nx; fb = brk_resid(b, &c);
        }
        if (!(fabs(fb) < 1e-8 && lo_ok < b && b < hi_ok)) continue;
        if (exact) {
            c.exact = 1;
            a = b; b = a + 1e-3;
            fa = brk_resid(a, &c); fb = brk_resid(b, &c);
            for (it = 0; it < 3; it++) {
                den = fb - fa;
                nx = b - (fabs(den) > 1e-300 ? fb * (b - a) / den : 0.0);
                a = b; fa = fb;
                b = nx; fb = brk_resid(b, &c);
            }
            if (!(fabs(fb) < 1e-8 && lo_ok < b && b < hi_ok)) continue;
        }
        lat[nbk] = b; left[nbk] = fl; right[nbk] = fr; nbk++;
    }
    return nbk;
}

typedef struct { double lat; int left, right; } Cut;

static int spot_rings(double theta_s, double rho, int nbk, const double *blat, const int *bl, const int *br,
                      double *theta, double *wth, double *dlam)
{
    Cut cuts[12], tmp;
    int nc = 0, b, c, k, nr = 0;
    double lo = theta_s - rho, hi = theta_s + rho, first, last;
    if (lo < 0.0) { cuts[nc++] = (Cut){0.0, 0, 0}; cuts[nc++] = (Cut){-lo, 0, 1}; }
    else cuts[nc++] = (Cut){lo, 1, 1};
    if (hi > M_PI) { cuts[nc++] = (Cut){TWO_PI - hi, 1, 0}; cuts[nc++] = (Cut){M_PI, 0, 0}; }
    else cuts[nc++] = (Cut){hi, 1, 1};
    first = cuts[0].lat; last = cuts[nc - 1].lat;
    for (b = 0; b < nbk; b++) {
        int clash = 0;
        if (!(first + 1e-6 < blat[b] && blat[b] < last - 1e-6)) continue;
        for (c = 0; c < nc; c++) if (fabs(blat[b] - cuts[c].lat) <= 1e-6) clash = 1;
        if (!clash) cuts[nc++] = (Cut){blat[b], bl[b], br[b]};
    }
    for (b = 1; b < nc; b++) {                      /* sort by latitude */
        tmp = cuts[b];
        for (c = b - 1; c >= 0 && cuts[c].lat > tmp.lat; c--) cuts[c + 1] = cuts[c];
        cuts[c + 1] = tmp;
    }
    for (c = 0; c + 1 < nc; c++) {
        double a = cuts[c].lat, e = cuts[c + 1].lat;
        int sa = cuts[c].right, sb = cuts[c + 1].left;
        if (e - a < 1e-9) continue;
        for (k = 0; k < NPP; k++) {
            double g = GLG[k], w = GLW[k], s, ds, t, cd;
            if (sa && sb) { s = sin(0.5 * M_PI * g); s *= s; ds = 0.5 * M_PI * sin(M_PI * g); }
            else if (sb) { s = 1.0 - (1.0 - g) * (1.0 - g); ds = 2.0 * (1.0 - g); }
            else if (sa) { s = g * g; ds = 2.0 * g; }
            else { s = g; ds = 1.0; }
            t = a + (e - a) * s;
            theta[nr] = t;
            wth[nr] = w * ds * (e - a);
            cd = (cos(rho) - cos(t) * cos(theta_s)) / fmax(sin(t) * sin(theta_s), 1e-300);
            dlam[nr] = acos(clip(cd, -1.0, 1.0));
            nr++;
        }
    }
    return nr;
}

/* ---------------------------------------------------------------- a ring's end directions */
typedef struct { const double *cth; double cg, sg, ci, ct, st, target; } EndCtx;

static double f_angle(double x, const EndCtx *e) { return mono(swept(e->cth, x)) - e->target; }

static double f_edge(double x, const EndCtx *e)       /* mu at direction x along the ring's own turn */
{
    double C = cos(swept(e->cth, x)), sp = sqrt(fmax(1.0 - C * C, 1e-300));
    return x * e->cg + sqrt(fmax(1.0 - x * x, 0.0)) * e->sg * (e->ci - e->ct * C) / (e->st * sp);
}

/* secant steps from two starts; if they do not meet the tolerance, a bracketed solve. Returns 0 if solved. */
static int solve(double (*f)(double, const EndCtx *), const EndCtx *e, double xa, double xb, double lo, double hi, double *root)
{
    double fa = f(xa, e), fb = f(xb, e), xn, flo, fhi;
    int it;
    for (it = 0; it < 12; it++) {
        if (fabs(fb) < 1e-13) { *root = xb; return 0; }
        if (fb == fa) break;
        xn = clip(xb - fb * (xb - xa) / (fb - fa), lo, hi);
        xa = xb; fa = fb;
        xb = xn; fb = f(xb, e);
    }
    if (fabs(fb) < 1e-13) { *root = xb; return 0; }
    flo = f(lo, e); fhi = f(hi, e);
    if (flo == 0.0) { *root = lo; return 0; }
    if (fhi == 0.0) { *root = hi; return 0; }
    if ((flo > 0.0) == (fhi > 0.0)) return 1;                 /* no change of sign: not bracketed */
    for (it = 0; it < 200 && hi - lo > 1e-16; it++) {
        double mid = 0.5 * (lo + hi), fm = f(mid, e);
        if (fm == 0.0) { *root = mid; return 0; }
        if ((fm > 0.0) == (flo > 0.0)) { lo = mid; flo = fm; } else { hi = mid; fhi = fm; }
    }
    *root = 0.5 * (lo + hi);
    return 0;
}

/* ---------------------------------------------------------------- the star's bin counts */
/* S[b * NE + e]: the integral over phase bin b (in turns) of the flux at energy e (photons / cm^2 / s / keV, times turns).
 * Returns 0, or a negative code: -1 outside the covered range, -2 a second path of light, -3 an end direction not found. */
static int spot_counts(double mass, double Req_km, double nu_hz, double D_cm, double incl, double theta_s, double rho, double kT,
                       const double *edges_turns, double phase0, double *S)
{
    double theta[MAXR], wth[MAXR], dlam[MAXR];        /* on this call's own stack, so that threads never share them */
    double GM = mass * XGM * C_CM * C_CM, L = 2.0 * mass * XGM, Req = Req_km * 1.0e5, Om = TWO_PI * nu_hz;
    double o2 = (Om * Om * Req * Req * Req / GM) * (-0.788 + 1.030 * GM / (Req * C_CM * C_CM));
    double u_eq = L / Req, u_pole = L / (Req * (1.0 + o2)), u_min = fmin(u_eq, u_pole), u_max = fmax(u_eq, u_pole);
    double h_max, ci = cos(incl), si = sin(incl), lam = -TWO_PI * phase0;
    double blat[4], edges[MAXB + 1], cth[MAXZ], cta[MAXZ], e_min = fmin(EN_E[0], EN_E[NE - 1]);
    int bl[4], br[4], nbk, nr, r, j, b, e, d, h, i, k;
    if (!(mass > 0.0 && Req_km > 0.0 && D_cm > 0.0 && kT > 0.0 && incl >= 0.0 && incl <= M_PI && rho > 0.0 && rho < M_PI
          && theta_s > 0.0 && theta_s < M_PI)) return -1;
    if (!(FIT_U0 < u_min && u_max < FIT_U1 && fabs(o2) < 0.5)) return -1;
    h_max = fabs(o2) / sqrt(1.0 - u_max) / (1.0 - fabs(o2));
    if (!(h_max < FIT_TILT)) return -1;
    nbk = visibility_breaks(L, Req, o2, incl, h_max > 0.05, blat, bl, br);
    if (8 * NPP > MAXR) return -1;                    /* at most 7 panels of NPP rings */
    nr = spot_rings(theta_s, rho, nbk, blat, bl, br, theta, wth, dlam);
    for (b = 0; b <= NB; b++) edges[b] = edges_turns[b] * TWO_PI;
    for (i = 0; i < NB * NE; i++) S[i] = 0.0;
    for (r = 0; r < nr; r++) {
        double ct = cos(theta[r]), st = sin(theta[r]), R = Req * (1.0 + o2 * ct * ct), u = L / R, f = sqrt(1.0 - u);
        double hh = -2.0 * Req * o2 * ct * st / (R * f), rt = sqrt(1.0 + hh * hh), cg = 1.0 / rt, sg = hh / rt, asg = fabs(sg);
        double beta = Om * R * st / (C_CM * f), gam = 1.0 / sqrt(1.0 - beta * beta), gi = 1.0 / gam;
        double A = ci * ct, B = si * st, trad = (Req - R + L * log((Req - L) / (R - L))) / C_CM;
        double area = R * R * st * rt * gam;
        double x_floor = -(asg + 0.02), x_crit = -sqrt(fmax(1.0 - 6.75 * u * u * (1.0 - u), 0.0));
        double c_hi, g0, g1, th_floor, g_floor, x_hi, x_far = 0.0, x_mu = 0.0, x_lo, mu_near, xm, xh, scale;
        double x[MAXN], th[MAXN], tau[MAXN], Cx[MAXN], phi[MAXN], sphi[MAXN], mu[MAXN], dly[MAXN];
        double delta[2][MAXN], wg[2][MAXN], a[2][MAXN], Np[2][MAXN][MAXE];
        double Cd[MAXD], phid[MAXD], Dd[MAXD], tabm[MAXD], tabp[MAXD], a_start;
        double Q[2][MAXB + 1], mq[2][MAXB + 1], lift[2][MAXB + 1], T[2][2][MAXB + 1][MAXN + 1], XI[2][2][MAXB + 1];
        int reach0, reach1, turning, full;
        char have[2][2][MAXB + 1];
        EndCtx ctx;
        if (x_floor < x_crit + 0.01) return -2;
        ring_series(u, cth, cta);
        /* No visible ray sweeps more than Theta(u, -|sin(tilt)|); a second path to some point of the ring needs a sweep of
         * at least 2 pi minus the angle of the ring's farthest point. If that can happen the star is refused. */
        if (swept(cth, -asg) >= TWO_PI - acos(clip(A - B, -1.0, 1.0)) - 1e-9) return -2;
        if (B < 1.0e-6) {                                   /* seen from along the axis: a steady flux from every longitude */
            double psi_c = acos(clip(A, -1.0, 1.0)), lo_x = x_floor, hi_x = 1.0, xc, xp, xmm, J, Dc, mu_c, dl_c, flux;
            int reach = swept(cth, lo_x) > psi_c;
            for (i = 0; i < 60; i++) { double mid = 0.5 * (lo_x + hi_x); if (swept(cth, mid) > psi_c) lo_x = mid; else hi_x = mid; }
            xc = 0.5 * (lo_x + hi_x);
            xp = fmin(xc + 1e-5, 1.0); xmm = xc - 1e-5;
            J = (xp - xmm) / fmax(fabs(cos(swept(cth, xp)) - cos(swept(cth, xmm))), 1e-300);
            Dc = ci * st / fmax(sin(psi_c), 1e-300);
            mu_c = xc * cg + sqrt(fmax(1.0 - xc * xc, 0.0)) * sg * clip(Dc, -1.0, 1.0);
            dl_c = gi;
            if (reach && mu_c > 0.0) {
                flux = wth[r] * area / (D_cm * D_cm) * 2.0 * dlam[r] * J * dl_c * dl_c * gi * mu_c;
                for (e = 0; e < NE; e++) {
                    double Ep = EN_E[e] / (f * dl_c), Ne = XBB * Ep * Ep / expm1(Ep / kT);
                    for (b = 0; b < NB; b++) S[b * NE + e] += (edges[b + 1] - edges[b]) / TWO_PI * flux * Ne;
                }
            }
            continue;
        }
        /* the two ends of the range of x, and the edge-on direction where part of the turn is hidden */
        c_hi = fmin(A + B, 1.0);
        g0 = 1.0 - c_hi;
        g1 = 1.0 - fmax(A - B, -1.0);
        th_floor = swept(cth, x_floor);
        g_floor = mono(th_floor);
        reach0 = g0 < g_floor - 1e-12;
        reach1 = g1 < g_floor - 1e-12;
        if (!reach0) continue;
        ctx.cth = cth; ctx.cg = cg; ctx.sg = sg; ctx.ci = ci; ctx.ct = ct; ctx.st = st;
        if (c_hi >= 1.0 - 1e-14) x_hi = 1.0;
        else {
            double xa = clip(1.0 - (1.0 - u) * g0, x_floor + 0.02, 1.0 - 1e-9);
            ctx.target = g0;
            if (solve(f_angle, &ctx, xa, clip(xa - 0.01, x_floor + 0.005, 1.0), x_floor, 1.0, &x_hi)) return -3;
        }
        turning = 0;
        if (reach1) {
            double xa = clip(1.0 - (1.0 - u) * g1, x_floor + 0.02, 1.0 - 1e-9), sigma = sin(theta[r] + incl) >= 0.0 ? 1.0 : -1.0;
            ctx.target = g1;
            if (solve(f_angle, &ctx, xa, clip(xa - 0.01, x_floor + 0.005, 1.0), x_floor, 1.0, &x_far)) return -3;
            turning = x_far * cg + sigma * sg * sqrt(fmax(1.0 - x_far * x_far, 0.0)) > 0.0;
        }
        mu_near = x_hi * cg + sqrt(fmax(1.0 - x_hi * x_hi, 0.0)) * sg * (theta[r] > incl ? 1.0 : (theta[r] < incl ? -1.0 : 0.0));
        x_mu = asg > 0.0 ? x_hi : 0.0;                      /* a round star's edge-on direction is x = 0 exactly */
        if (!turning && asg > 0.0 && mu_near > 0.0) {
            /* searched only between the directions of the farthest and the nearest point of the turn */
            double lo2 = fmax(-asg - 1e-9, reach1 ? x_far : x_floor), hi2 = fmin(asg + 1e-9, x_hi);
            if (hi2 > lo2 && solve(f_edge, &ctx, lo2, hi2, lo2, hi2, &x_mu)) return -3;
        }
        x_lo = turning ? x_far : fmin(x_mu, x_hi);
        if (!(mu_near > 0.0) || !(x_hi - x_lo > 1e-9)) continue;
        /* the points: x = xm + xh cos(t) */
        xm = 0.5 * (x_hi + x_lo); xh = 0.5 * (x_hi - x_lo);
        for (j = 0; j < NN; j++) {
            double cphi, spsi, salp, k1, bk;
            x[j] = xm + xh * CCOS[j];
            swept_delay(cth, cta, x[j], &th[j], &tau[j]);
            if (!(th[j] < M_PI + 1e-9)) return -2;        /* cannot happen after the test above; kept as a guard */
            Cx[j] = cos(th[j]);
            cphi = clip((Cx[j] - A) / B, -1.0, 1.0);
            sphi[j] = sqrt(fmax(1.0 - cphi * cphi, 1e-300));
            phi[j] = atan2(sphi[j], cphi);
            spsi = sin(th[j]);
            salp = sqrt(fmax(1.0 - x[j] * x[j], 0.0));
            mu[j] = fmax(x[j] * cg + salp * sg * (ci - ct * Cx[j]) / (st * spsi), 0.0);
            dly[j] = Om * (tau[j] * R / C_CM + trad);
            k1 = salp * si * sphi[j] / spsi;
            bk = beta * k1;
            delta[0][j] = gi / (1.0 - bk);                 /* before closest approach: the surface comes toward the observer */
            delta[1][j] = gi / (1.0 + bk);
            a[0][j] = dly[j] - phi[j];
            a[1][j] = dly[j] + phi[j];
            for (h = 0; h < 2; h++) {
                double dl = delta[h][j], s_ = 1.0 / (f * dl * kT);
                wg[h][j] = xh * CSIN[j] / (B * sphi[j]) * mu[j] * gi * dl * dl;
                if (e_min * s_ < 1.0e-3)              /* very low energy against the temperature: exp(y) - 1 would lose digits */
                    for (e = 0; e < NE; e++) { double y = EN_E[e] * s_; Np[h][j][e] = XBB * (y * kT) * (y * kT) / expm1(y); }
                else
                    for (e = 0; e < NE; e++) { double y = EN_E[e] * s_; Np[h][j][e] = XBB * (y * kT) * (y * kT) / (exp(y) - 1.0); }
            }
        }
        full = dlam[r] >= M_PI - 1e-12;
        scale = area * wth[r] / (D_cm * D_cm * TWO_PI);
        if (full) {
            for (e = 0; e < NE; e++) {
                double tot = 0.0;
                for (h = 0; h < 2; h++) for (j = 0; j < NN; j++) tot += CQ0[j] * wg[h][j] * Np[h][j][e];
                for (b = 0; b < NB; b++) S[b * NE + e] += scale * (edges[b + 1] - edges[b]) * tot;
            }
            continue;
        }
        /* arrival angle on the even grid */
        for (d = 0; d < ND; d++) { Cd[d] = 0.0; Dd[d] = 0.0; }
        for (j = 0; j < NN; j++) {
            const double *w = CTOD + (size_t)j * ND;
            double cj = Cx[j], dj = dly[j];
            for (d = 0; d < ND; d++) { Cd[d] += cj * w[d]; Dd[d] += dj * w[d]; }
        }
        for (d = 0; d < ND; d++) Cd[d] = clip((Cd[d] - A) / B, -1.0, 1.0);
        Cd[0] = clip((c_hi - A) / B, -1.0, 1.0);
        if (turning) Cd[ND - 1] = -1.0;
        phid[ND - 1] = 0.0;                                 /* (set in the loop below; this line only quiets a compiler warning) */
        for (d = 0; d < ND; d++) phid[d] = acos(Cd[d]);
        a_start = Dd[ND - 1] - phid[ND - 1];
        for (d = 0; d < ND; d++) { tabm[d] = Dd[ND - 1 - d] - phid[ND - 1 - d]; tabp[d] = Dd[d] + phid[d]; }
        memset(have, 0, sizeof have);
        /* queries: the bin edges shifted by the two ends of the stretch */
        for (k = 0; k < 2; k++) for (i = 0; i <= NB; i++) {
            double z = edges[i] + lam + (k == 0 ? dlam[r] : -dlam[r]), m = floor((z - a_start) / TWO_PI);
            mq[k][i] = m;
            Q[k][i] = z - TWO_PI * m;
            lift[k][i] = m * (a_start + TWO_PI) + M_PI * m * (m - 1.0) + m * (Q[k][i] - a_start);
        }
        for (h = 0; h < 2; h++) {
            const double *tab = h == 0 ? tabm : tabp;
            double t_first = tab[0], t_last = tab[ND - 1];
            for (k = 0; k < 2; k++) for (i = 0; i <= NB; i++) {
                double qf = Q[k][i], y0, y1, frac, pos;
                int lo = 0, hi = ND - 1, mid;
                if (qf <= t_first) { XI[h][k][i] = h == 0 ? 1.0 : -1.0; continue; }       /* before this half-turn's arrivals */
                if (qf >= t_last) { XI[h][k][i] = h == 0 ? -1.0 : 1.0; continue; }        /* after them */
                while (hi - lo > 1) { mid = (lo + hi) / 2; if (tab[mid] < qf) lo = mid; else hi = mid; }   /* tab[lo] < qf <= tab[hi] */
                y0 = tab[lo]; y1 = tab[lo + 1];
                frac = y1 > y0 ? (qf - y0) / (y1 - y0) : 0.0;
                pos = (lo + frac) / (ND - 1.0);
                XI[h][k][i] = clip(h == 0 ? 1.0 - 2.0 * pos : 2.0 * pos - 1.0, -1.0, 1.0);
            }
        }
        for (b = 0; b < NB; b++) {
            double DQ = Q[0][b + 1] - Q[0][b] - Q[1][b + 1] + Q[1][b];
            double Dl = lift[0][b + 1] - lift[0][b] - lift[1][b + 1] + lift[1][b];
            double Dm = mq[0][b + 1] - mq[0][b] - mq[1][b + 1] + mq[1][b];
            double acc[MAXE];
            for (e = 0; e < NE; e++) acc[e] = 0.0;
            for (h = 0; h < 2; h++) {
                double DT[MAXN + 1], DQT[MAXN + 1], c0 = h == 0 ? 1.0 : 0.0, c1 = h == 0 ? -1.0 : 1.0;
                double *T00 = T[h][0][b + 1], *T01 = T[h][0][b], *T10 = T[h][1][b + 1], *T11 = T[h][1][b];
                double Dp[MAXN], DQp[MAXN], x00 = XI[h][0][b + 1];
                if ((x00 == 1.0 || x00 == -1.0) && XI[h][0][b] == x00 && XI[h][1][b + 1] == x00 && XI[h][1][b] == x00) {
                    /* all four queries before, or all after, this half-turn's arrivals: the primitive is 0 or q0 at all four */
                    double pc = x00 == 1.0 ? 1.0 : 0.0;
                    for (j = 0; j < NN; j++) { Dp[j] = 0.0; DQp[j] = pc * DQ * CQ0[j]; }
                } else {
                    for (k = 0; k < 2; k++) for (i = b; i <= b + 1; i++) if (!have[h][k][i]) {      /* Chebyshev values, once per query */
                        double xi = XI[h][k][i], *Tq = T[h][k][i];
                        Tq[0] = 1.0; Tq[1] = xi;
                        for (j = 2; j <= NN; j++) Tq[j] = 2.0 * xi * Tq[j - 1] - Tq[j - 2];
                        have[h][k][i] = 1;
                    }
                    for (k = 0; k <= NN; k++) {
                        DT[k] = T00[k] - T01[k] - T10[k] + T11[k];
                        DQT[k] = Q[0][b + 1] * T00[k] - Q[0][b] * T01[k] - Q[1][b + 1] * T10[k] + Q[1][b] * T11[k];
                    }
                    for (j = 0; j < NN; j++) { Dp[j] = 0.0; DQp[j] = 0.0; }
                    for (k = 0; k <= NN; k++) {
                        const double *gk = CG + (size_t)k * NN;
                        double dt = DT[k], dq = DQT[k];
                        for (j = 0; j < NN; j++) { Dp[j] += dt * gk[j]; DQp[j] += dq * gk[j]; }
                    }
                }
                for (j = 0; j < NN; j++) {
                    double Mk = c1 * (DQp[j] - a[h][j] * Dp[j]) + CQ0[j] * (c0 * DQ + Dl - Dm * a[h][j]);
                    double wj = Mk * wg[h][j];
                    const double *np_ = Np[h][j];
                    for (e = 0; e < NE; e++) acc[e] += wj * np_[e];
                }
            }
            for (e = 0; e < NE; e++) S[b * NE + e] += scale * acc[e];
        }
    }
    return 0;
}

/* ---------------------------------------------------------------- the public likelihood's definition */
static double log_cond(double B, const double *d, const double *s, int n, double SCALE)
{
    double x = 0.0;
    int j;
    for (j = 0; j < n; j++) {
        double c = SCALE * (s[j] + B);
        if (c > 0.0) x += d[j] * log(c) - c;
        else if (fabs(c) < 1e-12 && fabs(d[j]) < 1e-12) continue;
        else return HUGE_NEG;
    }
    return x;
}

static double dlog(double B, const double *d, const double *s, int n, double T)
{
    double x = 0.0;
    int j;
    for (j = 0; j < n; j++) if (d[j] > 0.0) { if (s[j] + B > 0.0) x += d[j] / (s[j] + B); else return HUGE_POS; }
    return x - T;
}

/* the public code's search for the best background rate and the limits of its integral (spec 9.5, 9.6) */
static void channel_limits(const double *d, const double *s, int n, double T, double s0, double s1, double epsilon, double sigmas,
                           double *lower_, double *upper_, double *B_for_)
{
    double SCALE = T / n, ssum = 0.0, dsum = 0.0, av_STAR, av_DATA, B_min = 0.0, B, g, se = 0.0, std_est, lower, upper, B_for;
    int j;
    for (j = 0; j < n; j++) { ssum += s[j]; dsum += d[j]; }
    av_STAR = ssum / n; av_DATA = dsum / T;
    B = av_DATA - av_STAR;
    if (B <= B_min) {
        double min_counts = -1.0;
        for (j = 0; j < n; j++) if (fabs(s[j]) < 1e-12) { min_counts = -2.0; break; }
        if (min_counts == -2.0)
            for (j = 0; j < n; j++) if ((min_counts == -2.0 && d[j] > 0.0) || (0.0 < d[j] && d[j] < min_counts)) min_counts = d[j];
        if (min_counts != -1.0) {
            if (min_counts == -2.0) { B = 0.0; B_min = 0.0; }
            else { B = 0.01 * min_counts / SCALE; B_min = 0.1 * B; }
        } else B = B_min;
    }
    g = dlog(B_min, d, s, n, T);
    if (g <= 0.0) B = B_min;
    else {
        double B_low = B_min, B_high = av_DATA, y, x2, dB, std;
        int counter = 0;
        for (j = 0; j < n; j++) if (s[j] < 0.0 && av_DATA - s[j] > B_high) B_high = av_DATA - s[j];
        if (B < B_low || B >= B_high) B = 0.5 * (B_low + B_high);
        y = 0.0; x2 = 0.0;
        for (j = 0; j < n; j++) if (d[j] > 0.0) { double r = d[j] / (s[j] + B); y += r; x2 += r / (s[j] + B); }
        dB = (y - T) / x2; std = sqrt(1.0 / x2);
        while (!(fabs(dB) <= epsilon * std) && counter < 1000) {
            double B_new;
            g = y - T;
            if (g > 0.0) B_low = B; else if (g < 0.0) B_high = B; else break;
            B_new = B + dB;
            if (!(B_low < B_new && B_new < B_high)) B_new = 0.5 * (B_low + B_high);
            B = B_new;
            y = 0.0; x2 = 0.0;
            for (j = 0; j < n; j++) if (d[j] > 0.0) { double r = d[j] / (s[j] + B); y += r; x2 += r / (s[j] + B); }
            dB = (y - T) / x2; std = sqrt(1.0 / x2);
            counter++;
        }
    }
    for (j = 0; j < n; j++) if (d[j] > 0.0) se += d[j] / ((s[j] + B) * (s[j] + B));
    std_est = se > 0.0 ? sqrt(1.0 / se) : 1e90;
    lower = B - sigmas * std_est; upper = B + sigmas * std_est;
    if (lower < B_min) lower = B_min;
    B_for = B;
    if (lower < s0) {
        lower = s0;
        if (upper < s0 && s1 > 0.0) { upper = s1; B_for = s0; }
        else if (upper < s0) { upper = s0 + sigmas * std_est; B_for = s0; }
    }
    if (upper > s1 && s1 > 0.0) { upper = s1; if (lower > s1) { lower = s0; B_for = s1; } }
    if (B_for < lower) B_for = lower; else if (B_for > upper) B_for = upper;
    if (dlog(upper, d, s, n, T) >= 0.0) B_for = upper; else if (dlog(lower, d, s, n, T) <= 0.0) B_for = lower;
    *lower_ = lower; *upper_ = upper; *B_for_ = B_for;
}

/* the step-by-step form for one channel (used for channels with a zero model rate in a bin, or without counts) */
static double channel_reference(const double *d, const double *s, int n, double T, double s0, double s1, double fact)
{
    double SCALE = T / n, dsum = 0.0, ssum = 0.0, lower, upper, B_for, A, c2 = 0.0, g1 = 0.0, step0, lo, hi, step, val = 0.0, prev = -1.0;
    const double *g = getn("gl_g", 20)->v, *w = getn("gl_w", 20)->v;
    int j, panels, p, k;
    for (j = 0; j < n; j++) { dsum += d[j]; ssum += s[j]; }
    if (fabs(dsum / T) < 1e-12 && fabs(ssum / n) < 1e-12) {
        lower = fmax(0.0, s0); upper = 10.0 / T;
        if (upper > s1 && s1 > 0.0) upper = s1;
        return log((exp(-lower * T) - exp(-upper * T)) / T) + fact;
    }
    channel_limits(d, s, n, T, s0, s1, 1.0e-3, 10.0, &lower, &upper, &B_for);
    if (dsum == 0.0) return -SCALE * ssum - T * lower + log(-expm1(-T * (upper - lower)) / T) + fact;
    A = log_cond(B_for, d, s, n, SCALE);
    if (A < -1.0e299) return HUGE_NEG;
    for (j = 0; j < n; j++) if (d[j] > 0.0) { c2 += d[j] / ((s[j] + B_for) * (s[j] + B_for)); g1 += d[j] / (s[j] + B_for); }
    g1 = fabs(g1 - T);
    step0 = upper - lower;
    if (c2 > 0.0 && c2 < 1.0e299) step0 = fmin(step0, 1.0 / sqrt(c2));
    if (g1 > 0.0 && g1 < 1.0e299) step0 = fmin(step0, 1.0 / g1);
    lo = hi = B_for;
    step = step0;
    while (hi < upper) { hi = fmin(upper, hi + step); if (A - log_cond(hi, d, s, n, SCALE) >= 60.0) break; step *= 2.0; }
    step = step0;
    while (lo > lower) { lo = fmax(lower, lo - step); if (A - log_cond(lo, d, s, n, SCALE) >= 60.0) break; step *= 2.0; }
    for (panels = 16; panels <= 16384; panels *= 2) {
        double width = (hi - lo) / panels;
        val = 0.0;
        for (p = 0; p < panels; p++) for (k = 0; k < 20; k++) val += width * w[k] * exp(log_cond(lo + (p + g[k]) * width, d, s, n, SCALE) - A);
        if (prev >= 0.0 && fabs(val - prev) <= 1e-12 * fabs(val)) break;
        prev = val;
    }
    if (panels > 16384 || !(val > 0.0)) return HUGE_NEG;   /* did not settle, or zero: refused (the Python raises here) */
    return log(val) + A + fact;
}

static const double *GQ[6], *WQ[6];
static const int NQ[6] = {10, 14, 20, 28, 40, 64};
static const double QLIM[6] = {0.5, 2.0, 8.0, 18.0, 45.0, 1e300};

static double lfun(double B, const double *d, const double *s, int n, double SCALE, double T, double ssc)
{
    double x = 0.0;
    int j;
    for (j = 0; j < n; j++) x += d[j] * log(SCALE * (s[j] + B));
    return x - T * B - ssc;
}

/* star[r * NB + b]: the model's mean count rate in channel r and bin b. Returns the public code's log-likelihood.
 * PER_CHANNEL, if not NULL, gets each channel's value (the test of hard channels); it belongs to the caller. */
static double marginal_loglike(const double *star, double *PER_CHANNEL)
{
    const double DROP = 36.0, CUT_FALL = 32.0, T = T_EXP;
    double total = 0.0, SCALE = T / NB;
    int r, j, k, n = NB;
    for (r = 0; r < NCH; r++) {                      /* the public code's early exit: the model far above the data */
        double src = 0.0, dsum = 0.0;
        for (j = 0; j < n; j++) { src += star[r * n + j]; dsum += COUNTS[r * n + j]; }
        src *= SCALE;
        if (src - 20.0 * sqrt(fmax(src, 0.0)) - dsum > 0.0) return HUGE_NEG;
    }
    for (r = 0; r < NCH; r++) {
        const double *d = COUNTS + (size_t)r * n, *s = star + (size_t)r * n;
        double s0 = SUPPORT[2 * r], s1 = SUPPORT[2 * r + 1], dsum = 0.0, ssum = 0.0, smin = 1e300;
        double lower, upper, Bs, gs = 0.0, cs = 0.0, root, lo, hi, ssc, ls, fall_hi = 0.0, fall_lo = 0.0, span, fall, rmax = 0.0, rho, res, val = 0.0;
        double rinv[MAXB];
        int it, usual;
        for (j = 0; j < n; j++) { dsum += d[j]; ssum += s[j]; if (s[j] < smin) smin = s[j]; }
        usual = smin > 1e-12 && dsum > 0.0;
        if (!usual) { double cr = channel_reference(d, s, n, T, s0, s1, FACT[r]); if (cr < -1.0e299) return HUGE_NEG; if (PER_CHANNEL) PER_CHANNEL[r] = cr; total += cr; continue; }
        channel_limits(d, s, n, T, s0, s1, 1.0e-3, 10.0, &lower, &upper, &Bs);
        ssc = SCALE * ssum;
        ls = 0.0;
        for (j = 0; j < n; j++) {
            double q = s[j] + Bs, rj = d[j] / q;
            rinv[j] = 1.0 / q;
            if (rinv[j] > rmax) rmax = rinv[j];
            gs += rj; cs += rj / q;
            ls += d[j] * log(SCALE * q);
        }
        gs -= T;
        ls += -T * Bs - ssc;
        root = sqrt(gs * gs + 2.0 * DROP * cs);
        lo = fmax(lower, Bs - (gs <= 0.0 ? (root - gs) / cs : 2.0 * DROP / (root + gs)));
        hi = fmin(upper, Bs + (gs >= 0.0 ? (root + gs) / cs : 2.0 * DROP / (root - gs)));
        for (it = 0; it < 60; it++) {
            if (!(hi < upper && fall_hi < CUT_FALL)) break;
            fall_hi = ls - lfun(hi, d, s, n, SCALE, T, ssc);
            if (fall_hi < CUT_FALL) hi = fmin(upper, Bs + 1.5 * (hi - Bs));
        }
        for (it = 0; it < 60; it++) {
            if (!(lo > lower && fall_lo < CUT_FALL)) break;
            fall_lo = ls - lfun(lo, d, s, n, SCALE, T, ssc);
            if (fall_lo < CUT_FALL) lo = fmax(lower, Bs - 1.5 * (Bs - lo));
        }
        span = fmax(hi - Bs, Bs - lo);
        fall = fmax(fabs(gs) * span + 0.5 * cs * span * span, fmax(fall_hi, fall_lo));
        rho = span * rmax;
        if (dsum * pow(rho, 9) / (9.0 * (1.0 - fmin(rho, 0.5))) < 1e-10 && fall <= 2.0) {
            /* narrow limits: the logarithms as a power series about Bs, eight terms */
            double coef[8], pw[MAXB];
            const double *g = GQ[0], *w = WQ[0];
            for (j = 0; j < n; j++) pw[j] = d[j] * rinv[j];
            coef[0] = -T;
            for (j = 0; j < n; j++) coef[0] += pw[j];
            for (k = 2; k <= 8; k++) {
                double sum = 0.0;
                for (j = 0; j < n; j++) { pw[j] *= rinv[j]; sum += pw[j]; }
                coef[k - 1] = sum * ((k % 2 == 0) ? -1.0 : 1.0) / k;
            }
            for (k = 0; k < NQ[0]; k++) {
                double x = (lo - Bs) + (hi - lo) * g[k], acc = coef[7] * x;
                for (j = 6; j >= 0; j--) acc = (acc + coef[j]) * x;
                val += (hi - lo) * w[k] * exp(acc);
            }
        } else {
            int cls = 0;
            const double *g, *w;
            while (cls < 5 && fall > QLIM[cls]) cls++;
            g = GQ[cls]; w = WQ[cls];
            for (k = 0; k < NQ[cls]; k++) val += (hi - lo) * w[k] * exp(lfun(lo + (hi - lo) * g[k], d, s, n, SCALE, T, ssc) - ls);
        }
        if (!(val > 0.0)) return HUGE_NEG;                 /* an integral of zero: the public code's failure exit */
        res = log(val) + ls;
        if (PER_CHANNEL) PER_CHANNEL[r] = res + FACT[r];
        total += res + FACT[r];
    }
    return total;
}

/* ---------------------------------------------------------------- one whole working-out */
static double work_out(const double *p, double *S, double *star, int *code)
{
    double edges[MAXB + 1], incl = acos(p[3]), kT = XKB * pow(10.0, p[7]);
    int b, e, r;
    for (b = 0; b <= NB; b++) edges[b] = (double)b / NB;
    *code = spot_counts(p[0], p[1], FREQ, p[2] * XKPC, incl, p[5], p[6], kT, edges, 0.5 - p[4], S);
    if (*code) return HUGE_NEG;
    for (r = 0; r < NCH; r++) {
        const double *w = EN_W + (size_t)r * NE;
        for (b = 0; b < NB; b++) {
            double x = 0.0;
            for (e = 0; e < NE; e++) x += w[e] * S[b * NE + e];
            star[r * NB + b] = x;
        }
    }
    return marginal_loglike(star, NULL);
}

/* ---------------------------------------------------------------- many stars at once (also the shared library's three entry points) */
/* Once, before any batch: the tables. A missing or bad file ends the process with a message, as in the program. */
int nicer_load(const char *tables)
{
    load_tables(tables);
    return narr;
}

/* Once, or again for another setting, before a batch: rings per panel, points, energies (and the likelihood's quadratures). */
void nicer_choose(int npp, int n, int nE)
{
    int i;
    choose(npp, n, nE);
    for (i = 0; i < 6; i++) { GQ[i] = getn("gl_g", NQ[i])->v; WQ[i] = getn("gl_w", NQ[i])->v; }
}

/* The numbers each star's dump holds: the bin counts S (NB * NE), then the model's count rates (NCH * NB). */
int nicer_dump_size(void)
{
    return NB * NE + NCH * NB;
}

/* nstar stars: p holds eight numbers per star (as on the command line), lnl and code get one value per star (a refused star
 * gets its code and -1e300, as alone; its neighbours do not notice). bins, if not NULL, gets nicer_dump_size() numbers per star
 * (zeros for a refused star). The tables and the settings are only read here; each thread works in memory of its own. With
 * OpenMP (a compiler feature that shares a loop's turns among threads) the stars are handed out one at a time as threads
 * become free; without it the same loop runs on one thread. Returns 0, or 1 if working memory could not be had (those stars
 * get the code -9). */
int nicer_batch(int nstar, const double *p, int threads, double *lnl, int *code, double *bins)
{
    int i, failed = 0, finished = 0;
    size_t nS = (size_t)NB * NE, nR = (size_t)NCH * NB;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1)
#endif
    {
        double *S = (double *)malloc(nS * sizeof(double)), *star = (double *)malloc(nR * sizeof(double));
        int j, ok = S && star;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (j = 0; j < nstar; j++) {
            if (!ok) { lnl[j] = HUGE_NEG; code[j] = -9; continue; }
            lnl[j] = work_out(p + 8 * (size_t)j, S, star, &code[j]);
            if (bins) {
                double *o = bins + (size_t)j * (nS + nR);
                if (code[j]) memset(o, 0, (nS + nR) * sizeof(double));
                else { memcpy(o, S, nS * sizeof(double)); memcpy(o + nS, star, nR * sizeof(double)); }
            }
        }
        free(S); free(star);
        /* The OpenMP library's barrier at the end of this block already puts every thread's writes before what follows. The
         * race detector (make race) cannot see that barrier, as the library is not built for it; this counter shows it the
         * same order (each thread adds 1 after its writes, and the reading below waits for all of them). */
        __atomic_fetch_add(&finished, 1, __ATOMIC_RELEASE);
    }
    __atomic_load_n(&finished, __ATOMIC_ACQUIRE);
    for (i = 0; i < nstar; i++) if (code[i] == -9) failed = 1;
    (void)threads;
    return failed;
}

#ifndef NICER_LIBRARY
static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv)
{
    double *S = (double *)malloc(MAXB * MAXE * sizeof(double)), *star = (double *)malloc(NCHMAX * MAXB * sizeof(double));
    double p[8], v, best = 1e30, t0, tk = 0.0, tl = 0.0;
    int npp, n, nE, reps, i, k, code = 0, blocks = 7;
    if (!S || !star) { fprintf(stderr, "no memory\n"); return 2; }
    if (argc >= 8 && strcmp(argv[2], "batch") == 0) {
        /* nicer_fast <tables> batch <rings per panel> <points> <energies> <threads> <stars file> [dump file]: the file holds one
         * star per line (the eight numbers); one line per star is written: ln L to 17 significant digits (-inf for the public
         * code's early exit and for a refused star), then the code (0, or the refusal's code). The dump file, if named, gets
         * each star's bin counts and model rates, one star after another (zeros for a refused star). */
        FILE *fh = fopen(argv[7], "r");
        char line[4096];
        int cap = 1024, ns = 0, threads = atoi(argv[6]), *codes, nline = 0;
        double *ps = (double *)malloc((size_t)cap * 8 * sizeof(double)), *lnl, *bins = NULL;
        if (!fh) { fprintf(stderr, "cannot open %s\n", argv[7]); return 2; }
        if (threads < 1) { fprintf(stderr, "the number of threads must be 1 or more\n"); return 2; }
        while (fgets(line, sizeof line, fh)) {
            char *s = line, *end;
            int m;
            nline++;
            while (*s == ' ' || *s == '\t') s++;
            if (*s == '\n' || *s == '\r' || *s == 0) continue;
            if (ns == cap) { cap *= 2; ps = (double *)realloc(ps, (size_t)cap * 8 * sizeof(double)); }
            for (m = 0; m < 8; m++) {
                ps[8 * (size_t)ns + m] = strtod(s, &end);
                if (end == s) { fprintf(stderr, "%s, line %d: expected eight numbers\n", argv[7], nline); return 2; }
                s = end;
            }
            ns++;
        }
        fclose(fh);
#ifndef _OPENMP
        if (threads > 1) fprintf(stderr, "built without OpenMP: the stars are worked out on one thread\n");
#endif
        nicer_load(argv[1]);
        nicer_choose(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        lnl = (double *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(double)); codes = (int *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(int));
        if (argc > 8) bins = (double *)malloc((size_t)(ns > 0 ? ns : 1) * nicer_dump_size() * sizeof(double));
        if (!lnl || !codes || (argc > 8 && !bins)) { fprintf(stderr, "no memory\n"); return 2; }
        t0 = now_s();
        if (nicer_batch(ns, ps, threads, lnl, codes, bins)) { fprintf(stderr, "no memory for a thread's working memory\n"); return 2; }
        t0 = now_s() - t0;
        for (i = 0; i < ns; i++) { if (lnl[i] < -1.0e299) printf("-inf %d\n", codes[i]); else printf("%.17g %d\n", lnl[i], codes[i]); }
        if (bins) {
            fh = fopen(argv[8], "wb");
            fwrite(bins, sizeof(double), (size_t)ns * nicer_dump_size(), fh);
            fclose(fh);
        }
        fprintf(stderr, "%d stars on %d threads: %.6f s, %.1f stars per second\n", ns, threads, t0, ns / t0);
        return 0;
    }
    if (argc == 4 && strcmp(argv[2], "channels") == 0) {
        /* nicer_fast <tables> channels <file>: the file holds int32 channels, int32 bins, float64 exposure, then counts,
         * model rates, limits (2 per channel) and factorial sums; each channel's log-likelihood is written to <file>.out */
        FILE *fh = fopen(argv[3], "rb");
        int nch, nbin;
        double T, *cnt, *st, *sup, *fc, *outv, tot;
        char name[1024];
        if (!fh || fread(&nch, 4, 1, fh) != 1 || fread(&nbin, 4, 1, fh) != 1 || fread(&T, 8, 1, fh) != 1 || nbin > MAXB) { fprintf(stderr, "bad channel file\n"); return 2; }
        cnt = (double *)malloc((size_t)nch * nbin * 8); st = (double *)malloc((size_t)nch * nbin * 8);
        sup = (double *)malloc((size_t)nch * 16); fc = (double *)malloc((size_t)nch * 8); outv = (double *)calloc((size_t)nch, 8);
        if (fread(cnt, 8, (size_t)nch * nbin, fh) != (size_t)nch * nbin || fread(st, 8, (size_t)nch * nbin, fh) != (size_t)nch * nbin
            || fread(sup, 8, (size_t)nch * 2, fh) != (size_t)nch * 2 || fread(fc, 8, (size_t)nch, fh) != (size_t)nch) { fprintf(stderr, "short channel file\n"); return 2; }
        fclose(fh);
        load_tables(argv[1]);
        for (i = 0; i < 6; i++) { GQ[i] = getn("gl_g", NQ[i])->v; WQ[i] = getn("gl_w", NQ[i])->v; }
        NCH = nch; NB = nbin; T_EXP = T; COUNTS = cnt; SUPPORT = sup; FACT = fc;
        tot = marginal_loglike(st, outv);
        snprintf(name, sizeof name, "%s.out", argv[3]);
        fh = fopen(name, "wb");
        fwrite(outv, 8, (size_t)nch, fh);
        fclose(fh);
        if (tot < -1.0e299) printf("total = -inf\n"); else printf("total = %.9f\n", tot);
        return 0;
    }
    if (argc < 14) { fprintf(stderr, "usage: nicer_fast <tables> <rings per panel> <points> <energies> <repeats> mass radius distance cos_i phase_shift colat rho log10T [dump file]\n"); return 2; }
    load_tables(argv[1]);
    npp = atoi(argv[2]); n = atoi(argv[3]); nE = atoi(argv[4]); reps = atoi(argv[5]);
    for (i = 0; i < 8; i++) p[i] = atof(argv[6 + i]);
    nicer_choose(npp, n, nE);
    if (getenv("NICER_NO_UPPER_LIMIT")) {          /* the second public check's setting: no upper limit on the background */
        double *sup = (double *)malloc(2 * (size_t)NCH * sizeof(double));
        for (i = 0; i < NCH; i++) { sup[2 * i] = 0.0; sup[2 * i + 1] = -1.0; }
        SUPPORT = sup;
    }
    v = work_out(p, S, star, &code);
    if (code) { printf("refused: code %d\n", code); return 1; }
    if (v < -1.0e299) printf("ln L = -inf\n"); else printf("ln L = %.9f\n", v);
    printf("every bit of ln L: %a\n", v);              /* the same number, as a hexadecimal floating-point number */
    if (argc > 14) {
        FILE *fh = fopen(argv[14], "wb");
        fwrite(S, sizeof(double), (size_t)NB * NE, fh);
        fwrite(star, sizeof(double), (size_t)NCH * NB, fh);
        fclose(fh);
    }
    if (reps < 0) {                       /* a new star at every working-out: all eight numbers moved by up to 1 part in 1000 */
        double q[8], worst = 0.0;
        int nrep = -reps, bad = 0;
        for (k = 0; k < blocks; k++) {
            t0 = now_s();
            for (i = 0; i < nrep; i++) {
                int m;
                for (m = 0; m < 8; m++) q[m] = p[m] * (1.0 + 1.0e-3 * sin(1.7 * i + 2.3 * m + 0.4 * k));
                q[4] = 1.0e-3 * sin(0.9 * i + k);
                v = work_out(q, S, star, &code);
                if (code) bad++;
                if (v < worst) worst = v;
            }
            t0 = (now_s() - t0) / nrep;
            if (t0 < best) best = t0;
        }
        printf("time per working-out, a different star each time: %.4f ms (fastest of %d blocks of %d) ; refused %d ; lowest ln L %.3f\n",
               best * 1e3, blocks, nrep, bad, worst);
        return 0;
    }
    if (reps > 0) {
        double used = 0.0;
        for (k = 0; k < blocks; k++) {
            t0 = now_s();
            for (i = 0; i < reps; i++) { v = work_out(p, S, star, &code); used += v + S[i % (NB * NE)]; }
            t0 = (now_s() - t0) / reps;
            if (t0 < best) best = t0;
        }
        printf("sum of what the timed workings-out returned: %.6f\n", used);
        t0 = now_s();
        for (i = 0; i < reps; i++) { double edges[MAXB + 1]; int b; for (b = 0; b <= NB; b++) edges[b] = (double)b / NB;
            spot_counts(p[0], p[1], FREQ, p[2] * XKPC, acos(p[3]), p[5], p[6], XKB * pow(10.0, p[7]), edges, 0.5 - p[4], S); }
        tk = (now_s() - t0) / reps;
        t0 = now_s();
        for (i = 0; i < reps; i++) v = marginal_loglike(star, NULL);
        tl = (now_s() - t0) / reps;
        printf("time per working-out: %.4f ms (fastest of %d blocks of %d) ; the star's bin counts %.4f ms ; the likelihood step %.4f ms\n",
               best * 1e3, blocks, reps, tk * 1e3, tl * 1e3);
    }
    return 0;
}
#endif
