# -*- coding: utf-8 -*-
"""The fast way for a round or a flattened star (correctness version; see fast_kernel.py for the timed version).

J1  photons are counted at emission: the sum for one ring of latitude runs over x = cos(alpha), the direction the light leaves
    in, at a few Chebyshev points. Only the forward bending function is used.
J3  a stretch of the spot along a ring is added exactly, as the difference of the ring's cumulative count at its two ends.
J4  rings are placed by Gauss's rule across the spot.

THE TEST PROBLEMS' INPUTS (not derived by us; flagged in the plan):
  the flattened shape R(theta) = Req (1 + o2 cos^2 theta), o2 = Omega_bar^2 (-0.788 + 1.030 x), x = G M / (Req c^2),
  Omega_bar^2 = Omega^2 Req^3 / (G M)   [AlGendy and Morsink 2014; it was made to match general relativity's star models];
  the spot as a cap in the star's own angles; the black-body surface; the area convention of each test family.
OURS (from clocks and rulers; bend.py): the path, the delay, the clock-rate shift; special relativity for the surface's motion;
  the surface's tilt from the rulers: normal = (e_r - h e_theta)/sqrt(1 + h^2), h = R'/(R f), area R^2 sin(theta) sqrt(1 + h^2).
"""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bend import theta_tau_any
from consts import BB_CONST, C_CM_S, GM_SUN, PC_CM
from fast_sphere import _basis, bb_photons, gauss_rings

X_MIN = -0.45


class Star(object):
    def __init__(self, M_sun, Req_km, nu_hz, oblate=True):
        self.GM = GM_SUN * M_sun
        self.L = 2.0 * self.GM / C_CM_S**2
        self.Req = Req_km * 1.0e5
        self.Om = 2.0 * math.pi * nu_hz
        xc = self.GM / (self.Req * C_CM_S**2)
        ob2 = self.Om**2 * self.Req**3 / self.GM
        self.o2 = ob2 * (-0.788 + 1.030 * xc) if oblate else 0.0

    def radius(self, theta):
        c, s = math.cos(theta), math.sin(theta)
        return self.Req * (1.0 + self.o2 * c * c), -2.0 * self.Req * self.o2 * c * s

    def radial_delay(self, R):
        """Far-clock time for light to go radially from R to Req (zero for a round star): [F(Req) - F(R)]/c, F' = 1/f^2."""
        return (self.Req - R + self.L * math.log((self.Req - self.L) / (R - self.L))) / C_CM_S


class Ring(object):
    """One ring of latitude: cumulative count of the reference element (longitude 0) against its arrival angle."""

    def __init__(self, star, incl, theta, E_obs, kT, n=24, nd=769, ntab=260, beaming=None, gamma_area=True):
        self.empty = True
        R, dR = star.radius(theta)
        u = star.L / R
        f = math.sqrt(1.0 - u)
        h = dR / (R * f)
        cg, sg = 1.0 / math.sqrt(1.0 + h * h), h / math.sqrt(1.0 + h * h)
        st, ct_ = math.sin(theta), math.cos(theta)
        A, B = math.cos(incl) * ct_, math.sin(incl) * st
        if B < 1e-12:
            return
        beta = star.Om * R * st / (C_CM_S * f)
        gam = 1.0 / math.sqrt(1.0 - beta * beta)
        # the ring's own table of the forward function (the timed version replaces this by one fit for all stars)
        xt = np.concatenate([np.linspace(X_MIN, 0.05, 100, endpoint=False), np.linspace(0.05, 1.0, ntab - 100)])
        tht, taut = theta_tau_any(u, xt)
        ok = tht < math.pi - 1e-6
        xt, tht, taut = xt[ok], tht[ok], taut[ok]
        Ct = np.cos(tht)

        def mu_of(x, C):
            sp = np.sqrt(np.maximum(1.0 - C * C, 1e-300))
            cd = (math.cos(incl) - ct_ * C) / (st * sp)
            return x * cg + np.sqrt(np.maximum(1.0 - x * x, 0.0)) * sg * cd

        c_hi, c_far = A + B, A - B
        if c_hi <= Ct[0]:
            return
        x_hi = float(np.interp(c_hi, Ct, xt))
        if mu_of(x_hi, min(c_hi, 1.0 - 1e-15)) <= 0.0:
            return
        mu_t = mu_of(xt, Ct)
        turning_end = False
        if c_far > Ct[0]:
            x_far = float(np.interp(c_far, Ct, xt))
            below = (xt >= x_far) & (xt <= x_hi)
            if np.all(mu_t[below] > 0.0) and mu_of(x_far, c_far) > 0.0:
                x_lo, turning_end = x_far, True
        if not turning_end:
            idx = np.where((xt < x_hi) & (mu_t <= 0.0))[0]
            if idx.size == 0:
                x_lo = float(xt[0])
            else:
                k = idx[-1]                                   # last table point below x_hi where the surface is seen edge-on
                xa, xb = float(xt[k]), float(xt[k + 1])
                for _ in range(50):
                    xm = 0.5 * (xa + xb)
                    if mu_of(xm, float(np.interp(xm, xt, Ct))) > 0.0:
                        xb = xm
                    else:
                        xa = xm
                x_lo = xb
        xb_, hh = 0.5 * (x_hi + x_lo), 0.5 * (x_hi - x_lo)
        if hh < 5e-5:
            # A ring seen only within 1e-4 of its limb direction. Its two end directions come from table interpolation here
            # (errors near 1e-6), which is no longer small against the range, and a point can then fall outside the range
            # with an enormous weight. Its true contribution is below 1e-7 of a fully seen ring's, so it is left out.
            return
        self.empty = False
        tj, dct, td, to_dense, to_prim = _basis(n, nd)
        x = xb_ + hh * np.cos(tj)
        th, tau = theta_tau_any(u, x)
        Cx = np.cos(th)
        cphi = np.clip((Cx - A) / B, -1.0, 1.0)
        sphi = np.sqrt(np.maximum(1.0 - cphi * cphi, 1e-300))
        wgeo = hh * np.sin(tj) / (B * sphi)
        spsi = np.sin(th)
        salp = np.sqrt(np.maximum(1.0 - x * x, 0.0))
        mu = np.maximum(mu_of(x, Cx), 0.0)
        cphi_d = np.clip((to_dense @ Cx - A) / B, -1.0, 1.0)
        cphi_d[0] = min(1.0, (c_hi - A) / B)
        if turning_end:
            cphi_d[-1] = -1.0
        phi_d = np.arccos(cphi_d)
        delay_d = star.Om * ((to_dense @ tau) * R / C_CM_S + star.radial_delay(R))
        E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
        branches = []
        for s in (-1.0, 1.0):                                # s = -1: before closest approach (approaching side)
            cxi = -s * salp * math.sin(incl) * sphi / spsi
            delta = 1.0 / (gam * (1.0 - beta * cxi))
            Np = bb_photons(E_obs[None, :] / (f * delta[:, None]), kT)
            if beaming is not None:
                Np = Np * beaming(delta * mu)[:, None]
            w = (wgeo * delta**2 / gam * mu)[:, None] * Np
            branches.append((s * phi_d + delay_d, to_prim @ w))
        (am, cm), (ap, cp) = branches
        self.tot_minus = cm[-1]
        self.total = cm[-1] + cp[-1]
        self.a_start = am[-1]
        self.a_all = np.concatenate([am[::-1], ap[1:]])
        self.c_all = np.concatenate([(self.tot_minus[None, :] - cm)[::-1], self.tot_minus[None, :] + cp[1:]], axis=0)
        # area of the ring per unit longitude, in the convention of the test family
        self.area = R * R * st * math.sqrt(1.0 + h * h) * (gam if gamma_area else 1.0)

    def cumulative(self, q):
        q = np.asarray(q, dtype=float)
        two_pi = 2.0 * math.pi
        m = np.floor((q - self.a_start) / two_pi)
        qq = q - two_pi * m
        out = np.empty((q.size, self.c_all.shape[1]))
        for e in range(self.c_all.shape[1]):
            out[:, e] = np.interp(qq, self.a_all, self.c_all[:, e], right=self.total[e])
        return out + m[:, None] * self.total[None, :]


def light_curve(M_sun, Req_km, nu_hz, D_pc, incl, theta_s, rho, kT, E_obs, phases, n_rings=48, n=24, beaming=None,
                oblate=True, gamma_area=True):
    star = Star(M_sun, Req_km, nu_hz, oblate=oblate)
    D = D_pc * PC_CM
    theta, wth, dlam = gauss_rings(theta_s, rho, n_rings)
    a = np.asarray(phases) * 2.0 * math.pi
    E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
    out = np.zeros((a.size, E_obs.size))
    for th, w, dl in zip(theta, wth, dlam):
        ring = Ring(star, incl, th, E_obs, kT, n=n, beaming=beaming, gamma_area=gamma_area)
        if ring.empty:
            continue
        pref = ring.area * w / (D * D)
        if dl >= math.pi - 1e-12:
            out += pref * ring.total[None, :]
        else:
            out += pref * (ring.cumulative(a + dl) - ring.cumulative(a - dl))
    return out


SD1 = {"SD1a": (90.0, 0.01, 90.0, 1.0), "SD1b": (90.0, 1.0, 90.0, 1.0), "SD1c": (90.0, 0.01, 90.0, 200.0),
       "SD1d": (90.0, 1.0, 90.0, 200.0), "SD1e": (60.0, 1.0, 30.0, 400.0), "SD1f": (20.0, 1.0, 80.0, 400.0)}
OS1 = {"OS1a": (90.0, 0.01, 90.0, 600.0, None), "OS1b": (90.0, 1.0, 90.0, 600.0, None), "OS1c": (90.0, 0.01, 90.0, 200.0, None),
       "OS1d": (90.0, 1.0, 90.0, 1.0, None), "OS1e": (60.0, 1.0, 30.0, 600.0, None), "OS1f": (20.0, 1.0, 80.0, 600.0, None),
       "OS1g": (60.0, 1.0, 30.0, 600.0, "c"), "OS1h": (60.0, 1.0, 30.0, 600.0, "s"),
       "OS1i": (20.0, 1.0, 80.0, 600.0, "c"), "OS1j": (20.0, 1.0, 80.0, 600.0, "s")}
BEAM = {None: None, "c": lambda m: m * m, "s": lambda m: 1.0 - m * m}


def report(name, y, d):
    med = np.median(d[:, 1])
    e = (y - d[:, 1]) / med
    k = int(np.argmax(np.abs(e)))
    big = d[:, 1] > 0.01 * d[:, 1].max()
    print("%s  largest difference over the median %.2e at phase %.4f ; fractional, above 1%% of peak: %.2e ; r.m.s. over the median %.2e"
          % (name, abs(e[k]), d[k, 0], np.max(np.abs(y[big] / d[big, 1] - 1.0)), np.sqrt(np.mean(e**2))), flush=True)


if __name__ == "__main__":
    ref_dir = sys.argv[1]
    nr = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    for name, (ths, rho, inc, nu) in SD1.items():
        d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
        y = light_curve(1.4, 12.0, nu, 200.0, math.radians(inc), math.radians(ths), rho, 0.35, 1.0, d[:, 0], n_rings=nr,
                        oblate=False, gamma_area=False)[:, 0]
        report(name, y, d)
    for name, (ths, rho, inc, nu, bm) in OS1.items():
        d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
        y = light_curve(1.4, 12.0, nu, 200.0, math.radians(inc), math.radians(ths), rho, 0.35, 1.0, d[:, 0], n_rings=nr,
                        beaming=BEAM[bm], oblate=True, gamma_area=True)[:, 0]
        report(name, y, d)
