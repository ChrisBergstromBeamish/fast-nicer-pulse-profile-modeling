# -*- coding: utf-8 -*-
"""The timed kernel: photon counts of one hot spot in every phase bin and at every energy, for a round or a flattened star.

What is computed (derivation: NICER_SPEED_PLAN_v0.2.md section 3; Codex's reviews NRV1 T1 and NRV2 R3):
  for each ring of latitude r, each half-turn s (before or after closest approach), each Chebyshev point j in t, with
  x = cos(alpha) = xm + xh cos(t):
      w[r,s,j,E] = (xh sin t / (B sin phi)) (delta^2/gamma) mu N(E/(f delta); delta mu)      the emission-time count
      a[r,s,j]   = s phi + Omega (delay)                                                      the arrival angle
  A stretch of the spot (half-width dl) puts into the phase bin [b0, b1] the second difference
      K(b1+dl) - K(b0+dl) - K(b1-dl) + K(b0-dl),      K(z) = sum over emitted photons with a < z of (z - a),
  which is exact for the interpolated integrand: no patches, no curve to interpolate, no harmonic cut-off.
  K is linear in w, so the bins are filled by one matrix product at the end. The point where a = z is found by a table on an
  even grid and not refined: for the exact integrand an error in that point changes K only at second order (the summand
  vanishes there); for the interpolated integrand there is also a first-order term, the interpolation's residual at that
  point times the error (Codex, rounds 3 and 4). Measured by varying the grid alone from 65 to 2049 points: at the working
  129 points the bins change by 4e-9 of the mean bin and the quick example's ln L by 2e-6.
Light's path comes from clocks and rulers (see bend.py for the derivation). The flattened shape is the test problems' input.
"""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consts import BB_CONST, C_CM_S, GM_SUN

X_MIN = -0.45
ND = 129
_C0 = np.array([1.0, 0.0])[None, :, None, None]         # before / after closest approach
_C1 = np.array([-1.0, 1.0])[None, :, None, None]
_GLC = {}
_CH = {}
_TAB = {}


def _gl(n):
    if n not in _GLC:
        g, w = np.polynomial.legendre.leggauss(n)
        _GLC[n] = (0.5 * (g + 1.0), 0.5 * w)
    return _GLC[n]


def bend(u, x, n=40, need_delay=True):
    """Swept angle and delay (units R/c, behind a radial ray from the same radius) for arrays u = L/R and x = cos(alpha),
    x in (-1, 1]. Same integrals as bend.theta_tau_any; v = (xa/sqrt(a)) sinh(s) keeps the integrand smooth for every x."""
    u = np.asarray(u, dtype=float)
    x = np.asarray(x, dtype=float)
    shp = np.broadcast(u, x).shape
    u = np.broadcast_to(u, shp).ravel()
    x = np.broadcast_to(x, shp).ravel()
    xa = np.clip(np.abs(x), 1e-12, 1.0 - 1e-13)
    kb2 = (1.0 - xa * xa) / (1.0 - u)
    c1, c2, c3 = 2.0 - 3.0 * u, 1.0 - 3.0 * u, -u
    g, w = _gl(n)
    sq = np.sqrt(kb2 * c1)
    smax = np.arcsinh(sq / xa)
    s = smax[:, None] * g[None, :]
    es = np.exp(s)
    ei = 1.0 / es
    pref = (xa / sq)[:, None]
    v = pref * (0.5 * (es - ei))
    dv = pref * (0.5 * (es + ei)) * smax[:, None]
    w2 = v * v
    P = (xa * xa)[:, None] + kb2[:, None] * w2 * (c1[:, None] - c2[:, None] * w2 + c3[:, None] * w2 * w2)
    sP = np.sqrt(P)
    base = 2.0 * v * dv / sP
    th = np.sqrt(kb2) * (base @ w)
    ta = kb2 * ((base / (1.0 + sP)) @ w) if need_delay else None
    neg = x < 0.0
    if np.any(neg):
        xn, k2, un = xa[neg], kb2[neg], u[neg]
        d1, d2, d3 = c1[neg], c2[neg], c3[neg]
        st = xn * xn / (k2 * d1)
        for _ in range(6):
            F = k2 * st * (d1 + d2 * st + d3 * st * st) - xn * xn
            st = st - F / (k2 * (d1 + 2.0 * d2 * st + 3.0 * d3 * st * st))
        ss = st[:, None] * (1.0 - g * g)[None, :]
        gg = k2[:, None] * (d1[:, None] + d2[:, None] * (st[:, None] + ss) + d3[:, None] * (st[:, None] ** 2 + st[:, None] * ss + ss * ss))
        bs = 2.0 * np.sqrt(st)[:, None] / np.sqrt(gg)
        th[neg] += 2.0 * np.sqrt(k2) * (bs @ w)
        if need_delay:
            ta[neg] += 2.0 * ((bs / ((1.0 + ss) ** 2 * (1.0 - un[:, None] * (1.0 + ss)))) @ w)
    return th.reshape(shp), (ta.reshape(shp) if need_delay else None)


def _tables():
    """Three functions of u = L/R alone, tabulated once: the swept angle of a ray leaving along the surface (x = 0), its slope
    in x there, and cos of the swept angle at the lowest x used. They only place a cut in latitude and sort rings into kinds."""
    if not _TAB:
        ug = np.linspace(0.12, 0.62, 101)
        t0 = bend(ug, np.zeros_like(ug), n=96, need_delay=False)[0]
        tp = bend(ug, np.full_like(ug, 1e-4), n=96, need_delay=False)[0]
        tm = bend(ug, np.full_like(ug, -1e-4), n=96, need_delay=False)[0]
        tp2 = bend(ug, np.full_like(ug, 2e-3), n=96, need_delay=False)[0]
        tm2 = bend(ug, np.full_like(ug, -2e-3), n=96, need_delay=False)[0]
        with np.errstate(all="ignore"):
            tmin = bend(ug, np.full_like(ug, X_MIN), n=96, need_delay=False)[0]
        cm = np.cos(np.where(np.isfinite(tmin), np.minimum(tmin, math.pi), math.pi))   # no turning point: treat as out of range
        _TAB["t"] = (ug, t0, (tp - tm) / 2e-4, cm, (tp2 - 2.0 * t0 + tm2) / 4e-6)
        _TAB["lists"] = tuple(v.tolist() for v in _TAB["t"])
    return _TAB["t"]


def _tables_lists():
    """The same tables as plain Python lists (for the scalar arithmetic of the latitude cuts)."""
    _tables()
    return _TAB["lists"]


def _cheb(n):
    """Chebyshev tools on t in [0, pi], xi = 2 t/pi - 1, for n first-kind points."""
    if n not in _CH:
        k = np.arange(n)
        thj = math.pi * (k + 0.5) / n
        xi = np.cos(thj)
        tj = 0.5 * math.pi * (xi + 1.0)
        V = np.cos(np.outer(k, thj)) * (2.0 / n)             # coefficients c_k = sum_j V[k, j] f_j
        V[0] *= 0.5
        Bm = np.zeros((n + 1, n))                            # primitive from xi = -1
        for kk in range(n):
            if kk == 0:
                Bm[1, 0] += 1.0
            elif kk == 1:
                Bm[2, 1] += 0.25
                Bm[0, 1] += 0.25
            else:
                Bm[kk + 1, kk] += 0.5 / (kk + 1)
                Bm[kk - 1, kk] -= 0.5 / (kk - 1)
        sign = np.where(np.arange(n + 1) % 2 == 0, 1.0, -1.0)
        Bm[0, :] -= sign @ Bm
        G = 0.5 * math.pi * (Bm @ V)                         # p0_j(t) = sum_m T_m(xi) G[m, j]
        q = np.ones(n + 1) @ G                               # p0_j(pi)
        xid = np.linspace(-1.0, 1.0, ND)
        Td = np.cos(np.outer(np.arccos(xid), k))
        _CH[n] = (tj, np.cos(tj), np.sin(tj), G, q, xid, (Td @ V).T.copy(), np.arange(n + 1, dtype=float))
    return _CH[n]


# ---- one fit of light's path for all stars -----------------------------------------------------------------------------
# The swept angle and the delay depend on two numbers only, u = L/R and x = cos(alpha). With z = sqrt(1 - x) the functions
# Theta/z and tau/z^2 are smooth, so one two-variable Chebyshev fit holds them for every star in the range below. It is built
# once, from bend() itself, when first used (about a second). Measured against bend() with 96 points on 6000 random
# directions: largest error 9e-13 rad in the swept angle and 1.3e-12 R/c in the delay.
FIT_U = (0.15, 0.55)
FIT_XMIN = -0.12
FIT_NU, FIT_NZ = 32, 72
FIT_TILT = 0.08                       # the fit is used when every ring's tilt has |sin| below this (then x stays above FIT_XMIN)
_FIT = {}


def _path_fit():
    if not _FIT:
        (U0, U1), ZM = FIT_U, math.sqrt(1.0 - FIT_XMIN)
        ku, kz = np.arange(FIT_NU), np.arange(FIT_NZ)
        au, az = math.pi * (ku + 0.5) / FIT_NU, math.pi * (kz + 0.5) / FIT_NZ
        un = 0.5 * (U0 + U1) + 0.5 * (U1 - U0) * np.cos(au)
        zn = 0.5 * ZM * (1.0 + np.cos(az))
        th, ta = bend(un[:, None], (1.0 - zn * zn)[None, :], n=96)
        Vu = np.cos(np.outer(ku, au)) * (2.0 / FIT_NU)
        Vu[0] *= 0.5
        Vz = np.cos(np.outer(kz, az)) * (2.0 / FIT_NZ)
        Vz[0] *= 0.5
        _FIT["c"] = (Vu @ (th / zn[None, :]) @ Vz.T, Vu @ (ta / (zn * zn)[None, :]) @ Vz.T, ZM, ku.astype(float), kz.astype(float))
    return _FIT["c"]


def _ring_series(u):
    """Each ring's own one-variable series in z (coefficients of Theta/z and tau/z^2), from the two-variable fit."""
    Ct, Ca, ZM, ku, kz = _path_fit()
    Tu = np.cos(np.arccos(np.clip((2.0 * u - FIT_U[0] - FIT_U[1]) / (FIT_U[1] - FIT_U[0]), -1.0, 1.0))[:, None] * ku)
    return Tu @ Ct, Tu @ Ca, ZM, kz


def _series_eval(cth, cta, ZM, kz, x, need_delay=True):
    """Swept angle and delay at x (nr, m) from the rings' series."""
    z = np.sqrt(np.maximum(1.0 - x, 0.0))
    Tz = np.cos(np.arccos(np.clip(2.0 * z / ZM - 1.0, -1.0, 1.0))[..., None] * kz)      # (nr, m, FIT_NZ)
    th = z * np.matmul(Tz, cth[:, :, None])[..., 0]
    ta = z * z * np.matmul(Tz, cta[:, :, None])[..., 0] if need_delay else None
    return th, ta


def bb_photons(E, kT):
    return BB_CONST * E * E / np.expm1(E / kT)


def _panel(a, b, sa, sb, g, w):
    """Gauss nodes on [a, b]; an end marked True has square-root behaviour there and is smoothed by a squared variable."""
    if sa and sb:
        s = np.sin(0.5 * math.pi * g) ** 2
        ds = 0.5 * math.pi * np.sin(math.pi * g)
    elif sb:
        s = 1.0 - (1.0 - g) ** 2
        ds = 2.0 * (1.0 - g)
    elif sa:
        s = g ** 2
        ds = 2.0 * g
    else:
        s = g
        ds = np.ones_like(g)
    return a + (b - a) * s, w * ds * (b - a)


def spot_rings(theta_s, rho, n_per_panel, breaks=()):
    """Latitude nodes, weights and half-stretches for a circular spot. The range is cut where the spot starts to cover whole
    rings (a pole inside it) and at the given breaks (latitudes where part of the turn becomes hidden: the contribution of a
    ring has a corner there, with square-root behaviour on the hidden side)."""
    g, w = _gl(n_per_panel)
    lo, hi = theta_s - rho, theta_s + rho
    cuts = []                                                # (latitude, square-root on the left side, on the right side)
    if lo < 0.0:
        cuts += [(0.0, False, False), (-lo, False, True)]          # whole rings up to -lo, partial rings after it
    else:
        cuts += [(lo, True, True)]
    if hi > math.pi:
        cuts += [(2.0 * math.pi - hi, True, False), (math.pi, False, False)]
    else:
        cuts += [(hi, True, True)]
    first, last = cuts[0][0], cuts[-1][0]
    for b in breaks:
        if b is None:
            continue
        lat, left, right = (b, False, True) if isinstance(b, float) else b
        if first + 1e-6 < lat < last - 1e-6 and all(abs(lat - c[0]) > 1e-6 for c in cuts):
            cuts.append((lat, left, right))
    cuts.sort(key=lambda c: c[0])
    th, wt = [], []
    for (a, _, sa), (b, sb, _) in zip(cuts[:-1], cuts[1:]):
        if b - a < 1e-9:
            continue
        t_, w_ = _panel(a, b, sa, sb, g, w)
        th.append(t_)
        wt.append(w_)
    th = np.concatenate(th)
    wt = np.concatenate(wt)
    cd = (math.cos(rho) - np.cos(th) * math.cos(theta_s)) / np.maximum(np.sin(th) * math.sin(theta_s), 1e-300)
    return th, wt, np.arccos(np.clip(cd, -1.0, 1.0))


def _tab(u, tab):
    """Linear interpolation in a table on the even grid of u = L/R (plain Python numbers: this is called a few dozen times
    per working-out and must cost next to nothing)."""
    t = (u - 0.12) * 200.0                       # the grid is 0.12, 0.125, ..., 0.62
    i = min(max(int(t), 0), 99)
    return tab[i] + (t - i) * (tab[i + 1] - tab[i])


def visibility_breaks(L, Req, o2, incl, exact=None):
    """Latitudes at which a ring's visibility changes kind; the latitude range is cut there, because a ring's contribution has
    a corner there (with square-root or 3/2-power behaviour on one side).
    The far point of a ring's turn, at angle S = theta + i from the observer's direction (2 pi - S if S > pi), is seen edge-on
    where mu = x cos(tilt) + sigma sin(tilt) sqrt(1 - x^2) = 0, sigma = sign(sin S), i.e. at x = -sigma sin(tilt):
        Theta(u, -sin(tilt)) = theta + i            beyond it part of the turn is hidden
        Theta(u, +sin(tilt)) = 2 pi - theta - i     beyond it the whole turn is in view again (near the far pole)
    The near point, at angle |theta - i|, is seen edge-on at x = -sign(theta - i) sin(tilt):
        Theta(u, -sin(tilt)) = theta - i            beyond it the ring is never seen
        Theta(u, +sin(tilt)) = i - theta            before it the ring is never seen
    Each equation is solved by secant steps, with Theta near x = 0 from the tables to second order in the tilt (good to about
    2e-5 rad for a tilt of 0.03, to about 0.02 rad for a tilt of 0.2). Where the star's tilt is large and the caller gives
    the exact function, exact(u, x) for arrays, the contacts are then moved onto it by further secant steps.
    Returned as (latitude, corner side on the left, on the right) for spot_rings."""
    _, t0, t0x, _, t0xx = _tables_lists()
    two_pi = 2.0 * math.pi
    t_eq = _tab(L / Req, t0)
    kinds = ((-1.0, lambda T: T - incl, 0.0, math.pi - incl, (False, True)),
             (+1.0, lambda T: two_pi - incl - T, math.pi - incl, math.pi, (True, False)),
             (-1.0, lambda T: incl + T, incl, math.pi, (True, False)),
             (+1.0, lambda T: incl - T, 0.0, incl, (False, True)))

    def edge(th, sx):
        c, s = math.cos(th), math.sin(th)
        R = Req * (1.0 + o2 * c * c)
        u = L / R
        h = -2.0 * Req * o2 * c * s / (R * math.sqrt(1.0 - u))
        return u, sx * h / math.sqrt(1.0 + h * h)
    found = []
    for sx, lat_of, lo_ok, hi_ok, flags in kinds:
        lo_ok, hi_ok = max(lo_ok, 0.0) + 1e-4, min(hi_ok, math.pi) - 1e-4
        if hi_ok <= lo_ok:
            continue

        def F(th):
            u, x = edge(th, sx)
            return lat_of(_tab(u, t0) + _tab(u, t0x) * x + 0.5 * _tab(u, t0xx) * x * x) - th
        a = min(max(lat_of(t_eq), lo_ok), hi_ok)
        b = a + 0.02 if a + 0.02 < hi_ok else a - 0.02
        fa, fb = F(a), F(b)
        for _ in range(12):
            if abs(fb) < 1e-11 or fb == fa:
                break
            c = min(max(b - fb * (b - a) / (fb - fa), lo_ok), hi_ok)
            a, fa = b, fb
            b, fb = c, F(c)
        if abs(fb) < 1e-8 and lo_ok < b < hi_ok:
            found.append([b, sx, lat_of, lo_ok, hi_ok, flags])
    if exact is not None and found:
        sxs = np.array([f[1] for f in found])

        def Fx(th):
            ux = [edge(v, s_) for v, s_ in zip(th, sxs)]
            Th = exact(np.array([q[0] for q in ux]), np.array([q[1] for q in ux]))
            return np.array([f[2](T_) for f, T_ in zip(found, Th)]) - th
        a = np.array([f[0] for f in found])
        b = a + 1e-3
        fa, fb = Fx(a), Fx(b)
        for _ in range(3):
            den = fb - fa
            c = b - np.where(np.abs(den) > 1e-300, fb * (b - a) / np.where(den != 0.0, den, 1.0), 0.0)
            a, fa = b, fb
            b, fb = c, Fx(c)
        for f, v, r in zip(found, b, fb):
            f[0] = float(v) if (abs(r) < 1e-8 and f[3] < v < f[4]) else None
    return tuple((f[0],) + f[5] for f in found if f[0] is not None)


def _swept_exact(u, x, use_fit, n_gl=32):
    """Swept angle for arrays u, x of one shape, by the fit where it holds, otherwise by the integrals."""
    if use_fit and np.all(x >= FIT_XMIN):
        cth, cta, ZM, kz = _ring_series(u)
        return _series_eval(cth, cta, ZM, kz, x[:, None], need_delay=False)[0][:, 0]
    return bend(u, x, n=n_gl, need_delay=False)[0]


class TwoPaths(NotImplementedError):
    """The star is so compact (or so flattened) that light from behind reaches the observer by a second path around it.
    The program adds one path only and refuses such a star instead of returning a wrong number."""


ALIGN_B = 1.0e-6          # sin(inclination) sin(latitude) below which a ring is taken as seen from along the axis
END_TOL = 1.0e-12         # accepted residual of a ring's end directions
U_RANGE = (0.12, 0.62)    # L/R covered by the program


def _mono(th):
    """A smooth function of the swept angle that keeps growing past pi: 1 - cos up to pi, 3 + cos beyond."""
    c = np.cos(th)
    return np.where(th <= math.pi, 1.0 - c, 3.0 + c)


def spot_counts(M_sun, Req_km, nu_hz, D_cm, incl, theta_s, rho, kT, E_obs, edges, phase0=0.0, n_rings=10, n=16,
                oblate=True, gamma_area=True, beaming=None, rings=None, n_gl=32, debug=None, n_sec=3, fit=True):
    """Bin-integrated photon flux S[p, e]: the integral over phase bin p (in turns) of the flux at energy E_obs[e]
    (photons / cm^2 / s / keV, times turns). edges: bin edges in turns. phase0: the phase at which the spot's centre is closest.
    n_gl: Gauss points of the path integrals; n_sec: secant steps for each ring's end directions (their residuals are checked
    afterwards, and any that is not met is solved again inside a bracket). Measured on the public quick example's likelihood
    (ln L near -32552) and on its second, far more sensitive check star: 32 points change ln L by 7e-12 and 1e-6 against 96
    points; 3 steps by 7e-12 and 2e-6 against 8 steps.
    fit: use the one fit of light's path for all stars where it holds (see _path_fit); otherwise the integrals themselves.
    Covered: one path of light per surface point (TwoPaths is raised wherever a second path is possible), L/R inside U_RANGE."""
    GM = GM_SUN * M_sun
    L = 2.0 * GM / C_CM_S**2
    Req = Req_km * 1.0e5
    Om = 2.0 * math.pi * nu_hz
    o2 = (Om**2 * Req**3 / GM) * (-0.788 + 1.030 * GM / (Req * C_CM_S**2)) if oblate else 0.0
    u_eq, u_pole = L / Req, L / (Req * (1.0 + o2))
    if not (U_RANGE[0] < min(u_eq, u_pole) and max(u_eq, u_pole) < U_RANGE[1]):
        raise ValueError("L/R = %.3f to %.3f is outside the range the program covers (%.2f to %.2f)"
                         % (min(u_eq, u_pole), max(u_eq, u_pole), U_RANGE[0], U_RANGE[1]))
    h_max = abs(o2) / math.sqrt(1.0 - max(u_eq, u_pole)) / (1.0 - abs(o2))            # a bound on the tilt's tangent
    fit_all = fit and FIT_U[0] < min(u_eq, u_pole) and max(u_eq, u_pole) < FIT_U[1] and h_max < FIT_TILT
    if rings is None:
        # exact contacts only where the tilt is large (below a tilt of 0.05 the tables place the cuts to better than 1e-4 rad)
        exact = (lambda uu_, xx_: _swept_exact(uu_, xx_, fit_all, n_gl)) if h_max > 0.05 else None
        theta, wth, dlam = spot_rings(theta_s, rho, n_rings, breaks=visibility_breaks(L, Req, o2, incl, exact))
    else:
        theta, wth, dlam = rings
    E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
    edges = np.asarray(edges, dtype=float) * (2.0 * math.pi)
    nb = edges.size - 1
    ct_, st_ = np.cos(theta), np.sin(theta)
    R = Req * (1.0 + o2 * ct_ * ct_)
    u = L / R
    f = np.sqrt(1.0 - u)
    h = -2.0 * Req * o2 * ct_ * st_ / (R * f)
    rt = np.sqrt(1.0 + h * h)
    cg, sg = 1.0 / rt, h / rt
    beta = Om * R * st_ / (C_CM_S * f)
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    ci, si = math.cos(incl), math.sin(incl)
    A, B = ci * ct_, si * st_
    trad = (Req - R + L * np.log((Req - L) / (R - L))) / C_CM_S
    area = R * R * st_ * rt * (gam if gamma_area else 1.0)
    use_fit = fit_all and np.max(np.abs(sg)) < FIT_TILT
    asg = np.abs(sg)
    x_floor = -(asg + 0.02)                                        # no visible ray starts further inward than the tilt
    x_crit = -np.sqrt(np.maximum(1.0 - 6.75 * u * u * (1.0 - u), 0.0))   # inward of this a ray does not come back out
    if np.any(x_floor < x_crit + 0.01):
        raise TwoPaths("a ring's tilt (sin %.3f) reaches the directions that do not come back out (L/R %.3f)"
                       % (float(asg.max()), float(u.max())))
    if use_fit:
        cth, cta, ZM, kz = _ring_series(u)

    def swept(X, idx=None, delay=False):
        """Swept angle (and delay) at X of shape (k, rings) or, with idx, at X[m] for the rings idx[m]."""
        if idx is None:
            if use_fit:
                th, ta = _series_eval(cth, cta, ZM, kz, X.T, need_delay=delay)
                return th.T, (ta.T if delay else None)
            return bend(np.broadcast_to(u, X.shape), X, n=n_gl, need_delay=delay)
        if use_fit:
            th, ta = _series_eval(cth[idx], cta[idx], ZM, kz, X[:, None], need_delay=delay)
            return th[:, 0], (ta[:, 0] if delay else None)
        return bend(u[idx], X, n=n_gl, need_delay=delay)

    def tilt_term(X, th, idx=None):
        """mu at direction X with swept angle th (below pi): x cos(tilt) + sin(alpha) sin(tilt) cos(delta)."""
        c_g, s_g, c_t, s_t = (cg, sg, ct_, st_) if idx is None else (cg[idx], sg[idx], ct_[idx], st_[idx])
        C = np.cos(th)
        sp = np.sqrt(np.maximum(1.0 - C * C, 1e-300))
        return X * c_g + np.sqrt(np.maximum(1.0 - X * X, 0.0)) * s_g * (ci - c_t * C) / (s_t * sp)

    # ---- one path of light only (Codex's round-4 counterexample: a round star with L/R = 0.6 seen from its axis, a spot near
    # the far pole; the tangent ray sweeps more than pi, so a point at angle psi is also reached by a ray sweeping 2 pi - psi).
    # A ray that leaves the surface has x > -|sin(tilt)|, and the swept angle grows as x falls, so no visible ray sweeps more
    # than Theta(u, -|sin(tilt)|). The farthest point of a ring's turn is at angle psi_far; a second path to some point of the
    # ring needs a sweep of at least 2 pi - psi_far. If the largest sweep is below that, there is one path only; if not, the
    # program refuses the star.
    th_edge = swept(-asg[None, :])[0][0]
    psi_far = np.arccos(np.clip(A - B, -1.0, 1.0))
    if np.any(th_edge >= 2.0 * math.pi - psi_far - 1e-9):
        k_ = int(np.argmax(th_edge + psi_far))
        raise TwoPaths("light from behind the star reaches the observer by a second path: at latitude %.3f rad a visible ray can "
                       "sweep %.3f rad, and the far point of that ring is %.3f rad from the observer's direction" % (theta[k_], th_edge[k_], psi_far[k_]))
    # ---- rings seen from along the spin axis: every longitude gives the same steady flux
    S_al = np.zeros((nb, E_obs.size))
    aligned = B < ALIGN_B
    if np.any(aligned):
        ia = np.where(aligned)[0]
        psi_c = np.arccos(np.clip(A[ia], -1.0, 1.0))
        lo_x, hi_x = x_floor[ia].copy(), np.ones(ia.size)
        reach = swept(lo_x, ia)[0] > psi_c
        for _ in range(60):
            mid = 0.5 * (lo_x + hi_x)
            big = swept(mid, ia)[0] > psi_c
            lo_x, hi_x = np.where(big, mid, lo_x), np.where(big, hi_x, mid)
        xc = 0.5 * (lo_x + hi_x)
        dx = 1e-5
        xp, xm_ = np.minimum(xc + dx, 1.0), xc - dx
        J = (xp - xm_) / np.maximum(np.abs(np.cos(swept(xp, ia)[0]) - np.cos(swept(xm_, ia)[0])), 1e-300)
        Dc = ci * st_[ia] / np.maximum(np.sin(psi_c), 1e-300)
        mu_c = xc * cg[ia] + np.sqrt(np.maximum(1.0 - xc * xc, 0.0)) * sg[ia] * np.clip(Dc, -1.0, 1.0)
        dl_c = 1.0 / gam[ia]
        Nc = bb_photons(E_obs[None, :] / (f[ia] * dl_c)[:, None], kT)
        if beaming is not None:
            Nc = Nc * beaming(dl_c * mu_c)[:, None]
        flux = (wth[ia] * area[ia] / (D_cm * D_cm) * 2.0 * dlam[ia] * J * dl_c * dl_c / gam[ia] * np.maximum(mu_c, 0.0)
                * (reach & (mu_c > 0.0)))[:, None] * Nc
        S_al = ((edges[1:] - edges[:-1]) / (2.0 * math.pi))[:, None] * flux.sum(axis=0)[None, :]

    # ---- the ends of each ring's range of x. First the directions for the nearest and the farthest point of the turn (two
    # unknowns per ring, together). Then, on rings with a hidden part, the edge-on direction, searched only between those two
    # (outside them the angle to the observer is not on the ring's turn and mu has no meaning).
    def solve_rows(fun, Xa, Xb, lo, hi, need):
        """Secant steps from two starts for an array of unknowns; every needed one is checked, and those not met are solved
        again inside their bracket."""
        Fa, Fb = fun(Xa), fun(Xb)
        for it in range(n_sec + 6):
            if it >= n_sec and not np.any(need & ~(np.abs(Fb) < END_TOL)):
                break
            den = Fb - Fa
            Xn = np.clip(Xb - np.where(np.abs(den) > 1e-300, Fb * (Xb - Xa) / np.where(den != 0.0, den, 1.0), 0.0), lo, hi)
            Xa, Fa = Xb, Fb
            Xb = Xn
            Fb = fun(Xb)
        bad = need & ~(np.abs(Fb) < END_TOL)
        if np.any(bad):
            idx = np.where(bad)
            lo_x, hi_x = lo[idx].copy(), hi[idx].copy()
            f_lo, f_hi = fun(lo_x, idx), fun(hi_x, idx)
            if np.any((np.sign(f_lo) == np.sign(f_hi)) & (f_lo != 0.0) & (f_hi != 0.0)):
                raise ArithmeticError("an end direction of a ring could not be bracketed (latitudes %s rad)"
                                      % np.round(theta[idx[-1]], 4).tolist())
            for _ in range(60):
                mid = 0.5 * (lo_x + hi_x)
                same = np.sign(fun(mid, idx)) == np.sign(f_lo)
                lo_x, hi_x = np.where(same, mid, lo_x), np.where(same, hi_x, mid)
            Xb[idx] = 0.5 * (lo_x + hi_x)
        return Xb

    c_hi = np.minimum(A + B, 1.0)
    g_tgt = np.stack([1.0 - c_hi, 1.0 - np.maximum(A - B, -1.0)])              # 1 - cos of the nearest and farthest angle
    th_floor = swept(x_floor[None, :])[0][0]
    g_floor = _mono(th_floor)
    reach = g_tgt < g_floor[None, :] - 1e-12                                    # the angle can be reached at all
    g_use = np.where(reach, g_tgt, 0.5 * g_floor[None, :])                      # a harmless target where it cannot

    def f_angle(X, idx=None):
        if idx is None:
            return _mono(swept(X)[0]) - g_use
        return _mono(swept(X, idx[1])[0]) - g_use[idx]
    lo1 = np.stack([x_floor, x_floor])
    hi1 = np.ones((2, u.size))
    X1 = np.clip(1.0 - (1.0 - u) * g_use, lo1 + 0.02, 1.0 - 1e-9)
    X1 = solve_rows(f_angle, X1, np.clip(X1 - 0.01, lo1 + 0.005, hi1), lo1, hi1, reach)
    x_hi = np.where(c_hi >= 1.0 - 1e-14, 1.0, X1[0])
    x_far = X1[1]
    mu_near = x_hi * cg + np.sqrt(np.maximum(1.0 - x_hi * x_hi, 0.0)) * sg * np.sign(theta - incl)
    sigma = np.where(np.sin(theta + incl) >= 0.0, 1.0, -1.0)
    # the far point of the turn in view: then the whole turn is taken as in view and no edge-on direction is needed
    turning = reach[1] & (x_far * cg + sigma * sg * np.sqrt(np.maximum(1.0 - x_far * x_far, 0.0)) > 0.0)
    lo2 = np.maximum(-asg - 1e-9, np.where(reach[1], x_far, x_floor))
    hi2 = np.minimum(asg + 1e-9, x_hi)
    need2 = (asg > 0.0) & reach[0] & ~turning & (mu_near > 0.0) & ~aligned & (hi2 > lo2)
    x_mu = np.where(asg > 0.0, x_hi, 0.0)                          # a round star's edge-on direction is x = 0 exactly
    if np.any(need2):
        lo2, hi2 = np.where(need2, lo2, 0.0), np.where(need2, hi2, 0.0)

        def f_edge(X, idx=None):
            if idx is None:
                return tilt_term(X, swept(X[None, :])[0][0])
            return tilt_term(X, swept(X, idx[0])[0], idx[0])
        x_mu = np.where(need2, solve_rows(f_edge, lo2.copy(), hi2.copy(), lo2, hi2, need2), x_mu)
    x_lo = np.where(turning, x_far, np.minimum(x_mu, x_hi))
    keep = reach[0] & (mu_near > 0.0) & ~aligned & (x_hi - x_lo > 1e-9)
    if not np.any(keep):
        return S_al
    if not np.all(keep):
        sel = np.where(keep)[0]
        (theta, wth, dlam, ct_, st_, R, u, f, rt, cg, sg, beta, gam, A, B, trad, x_hi, x_lo, turning, c_hi, area) = [
            v[sel] for v in (theta, wth, dlam, ct_, st_, R, u, f, rt, cg, sg, beta, gam, A, B, trad, x_hi, x_lo, turning, c_hi, area)]
        if use_fit:
            cth, cta = cth[sel], cta[sel]
    nr = theta.size
    if debug is not None:
        debug.update(theta=theta, x_hi=x_hi, x_lo=x_lo, turning=turning, c_hi=c_hi, A=A, B=B, u=u)
    tj, costj, sintj, G, q0, xid, to_d, kk = _cheb(n)
    xm, xh = 0.5 * (x_hi + x_lo), 0.5 * (x_hi - x_lo)
    x = xm[:, None] + xh[:, None] * costj[None, :]                           # (nr, n)
    th, tau = _series_eval(cth, cta, ZM, kz, x) if use_fit else bend(u[:, None], x, n=n_gl)
    if not np.all(th < math.pi + 1e-9):                          # cannot happen after the test above; kept as a guard
        raise TwoPaths("a swept angle above pi on the first path (%.4f rad)" % float(np.nanmax(th)))
    Cx = np.cos(th)
    cphi = np.clip((Cx - A[:, None]) / B[:, None], -1.0, 1.0)
    sphi = np.sqrt(np.maximum(1.0 - cphi * cphi, 1e-300))
    phi = np.arctan2(sphi, cphi)
    spsi = np.sin(th)
    salp = np.sqrt(np.maximum(1.0 - x * x, 0.0))
    mu = np.maximum(x * cg[:, None] + salp * sg[:, None] * (ci - ct_[:, None] * Cx) / (st_[:, None] * spsi), 0.0)
    dly = Om * (tau * R[:, None] / C_CM_S + trad[:, None])                    # (nr, n)
    k1 = salp * si * sphi / spsi                                              # |cos xi|
    bk = beta[:, None] * k1
    gi = 1.0 / gam[:, None]
    delta = np.stack([gi / (1.0 + bk), gi / (1.0 - bk)], axis=1)              # s = -1: cos xi = +k1, approaching
    delta = delta[:, ::-1, :]                                                 # index 0: s = -1 (before closest approach)
    wg = (xh[:, None] * sintj[None, :] / (B[:, None] * sphi) * mu * gi)[:, None, :] * delta * delta
    if beaming is not None:
        wg = wg * beaming(delta * mu[:, None, :])
    Np = bb_photons(E_obs[None, None, None, :] / (f[:, None, None, None] * delta[..., None]), kT)
    a = np.stack([dly - phi, dly + phi], axis=1)                              # (nr, 2, n) arrival angle at the nodes
    if debug is not None:
        debug.update(x=x, Cx=Cx, cphi=cphi, sphi=sphi, mu=mu, wg=wg, Np=Np, a=a, delta=delta)
    # ---- arrival angle on the even grid in xi (first guess is all that is needed)
    Cd = np.clip((Cx @ to_d - A[:, None]) / B[:, None], -1.0, 1.0)            # (nr, ND); xi = -1 is closest approach
    Cd[:, 0] = np.clip((c_hi - A) / B, -1.0, 1.0)
    Cd[:, -1] = np.where(turning, -1.0, Cd[:, -1])
    phid = np.arccos(Cd)
    Dd = dly @ to_d
    a_start = Dd[:, -1] - phid[:, -1]
    a_mid = Dd[:, 0]
    a_end = Dd[:, -1] + phid[:, -1]
    two_pi = 2.0 * math.pi
    # ---- queries: bin edges shifted by the two ends of the stretch
    lam = -two_pi * phase0
    z = edges[None, None, :] + lam + np.array([1.0, -1.0])[None, :, None] * dlam[:, None, None]      # (nr, 2 ends, nb + 1)
    mq = np.floor((z - a_start[:, None, None]) / two_pi).reshape(nr, -1)
    Q = z.reshape(nr, -1) - two_pi * mq                                       # in [a_start, a_start + 2 pi)
    nq = Q.shape[1]
    # inverse interpolation for both half-turns at once: rows are made increasing and separated by an offset
    off = (np.arange(nr) * 100.0)[:, None]
    inc_minus = (Dd - phid)[:, ::-1]                                          # increasing in the reversed grid
    inc_plus = Dd + phid
    xis = np.empty((nr, 2, nq))
    for si_, (tab, rev) in enumerate(((inc_minus, True), (inc_plus, False))):
        flat = (tab + off).ravel()
        qf = (np.clip(Q, tab[:, :1], tab[:, -1:]) + off).ravel()
        idx = np.clip(np.searchsorted(flat, qf) - 1, 0, flat.size - 2)
        row = np.repeat(np.arange(nr), nq)
        loc = np.clip(idx - row * ND, 0, ND - 2)
        base = row * ND + loc
        y0, y1 = flat[base], flat[base + 1]
        frac = np.where(y1 > y0, (qf - y0) / np.where(y1 > y0, y1 - y0, 1.0), 0.0)
        pos = (loc + frac).reshape(nr, nq) / (ND - 1.0)                       # in [0, 1] along the increasing table
        xis[:, si_] = (1.0 - 2.0 * pos) if rev else (2.0 * pos - 1.0)
    xc = np.clip(xis, -1.0, 1.0)                                              # Chebyshev values by the three-term rule
    Tk1 = np.empty(xc.shape + (n + 1,))                                       # (nr, 2, nq, n + 1)
    Tk1[..., 0] = 1.0
    Tk1[..., 1] = xc
    x2 = 2.0 * xc
    for m_ in range(2, n + 1):
        Tk1[..., m_] = x2 * Tk1[..., m_ - 1] - Tk1[..., m_ - 2]
    # With p0 = Tk1 @ G (the integral of each node's cardinal function from t = 0 to the query's t*; a query outside the
    # half-turn's range has t* at the nearer end, so p0 is 0 or q0 there), the kernel of K at one query is
    #     (Q - a) part + q0 (lift - m a),      part = p0 after closest approach, q0 - p0 before it.
    # Everything is linear in the query's quantities, so the second difference D over the bin's two edges and the stretch's
    # two ends is taken first, on small arrays, and G is applied to the differences:
    #     D[kernel] = c1 (D[Q p0] - a D[p0]) + q0 (c0 D[Q] + D[lift] - a D[m]),   (c0, c1) = (1, -1) before, (0, 1) after.
    T5 = Tk1.reshape(nr, 2, 2, nb + 1, n + 1)                                 # (ring, half-turn, stretch end, edge, order)
    QT = Q.reshape(nr, 1, 2, nb + 1, 1) * T5
    DT = (T5[:, :, 0, 1:] - T5[:, :, 0, :-1]) - (T5[:, :, 1, 1:] - T5[:, :, 1, :-1])
    DQT = (QT[:, :, 0, 1:] - QT[:, :, 0, :-1]) - (QT[:, :, 1, 1:] - QT[:, :, 1, :-1])
    Dp = DT @ G                                                               # (nr, 2, nb, n)
    DQp = DQT @ G

    def D2(X):
        X = X.reshape(nr, 2, nb + 1)
        return ((X[:, 0, 1:] - X[:, 0, :-1]) - (X[:, 1, 1:] - X[:, 1, :-1]))[:, None, :, None]
    lift = mq * (a_start[:, None] + two_pi) + math.pi * mq * (mq - 1.0) + mq * (Q - a_start[:, None])
    an = a[:, :, None, :]
    Mk = _C1 * (DQp - an * Dp) + q0 * ((_C0 * D2(Q) + D2(lift)) - D2(mq) * an)
    full = dlam >= math.pi - 1e-12
    if np.any(full):
        Mk[full] = q0[None, None, None, :] * (edges[1:] - edges[:-1])[None, None, :, None]
    area = R * R * st_ * rt * (gam if gamma_area else 1.0)
    Mk *= (area * wth / (D_cm * D_cm * two_pi))[:, None, None, None]
    Mw = (Mk * wg[:, :, None, :]).transpose(2, 0, 1, 3).reshape(nb, -1)       # (nb, nr * 2 * n)
    return S_al + Mw @ Np.reshape(-1, E_obs.size)
