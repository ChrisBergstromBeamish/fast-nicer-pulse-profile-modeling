# -*- coding: utf-8 -*-
"""The fast way for a round star, first version (correctness before speed).

What it does differently from the usual way (NICER_SPEED_PLAN_v0.2.md, J1, J3, J4):
  J1  photons are counted at emission: for one ring of latitude the sum runs over cos(alpha), the direction the light leaves in,
      at a few Chebyshev points. Only the forward bending function is used: no direction search, no lensing factor.
      Emission-time weight (derived twice, by Claude and by Codex):  (delta^2 / gamma) mu N(E/(f delta); delta mu) per unit
      emission phase and unit longitude, times  d cos(alpha) / sqrt(B^2 - (C - A)^2).
  J3  a stretch of the spot along the ring is added exactly: the flux at arrival phase a is the difference of the ring's
      cumulative count at the two ends of the stretch (no patches, no interpolation of a curve, no cut-off).
  J4  rings are placed by Gauss's rule across the spot, with the spot's outline giving each ring's stretch.
"""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bend import theta_tau
from consts import BB_CONST, C_CM_S, GM_SUN, PC_CM

_BASIS = {}


def _basis(n, nd):
    """Cosine-series tools for n midpoint nodes t_j = pi (j + 1/2)/n on (0, pi), and a dense grid of nd points on [0, pi]."""
    key = (n, nd)
    if key not in _BASIS:
        tj = math.pi * (np.arange(n) + 0.5) / n
        k = np.arange(n)
        dct = np.cos(np.outer(k, tj)) * (2.0 / n)          # coefficients c_k = sum_j dct[k, j] f_j
        dct[0] *= 0.5
        td = np.linspace(0.0, math.pi, nd)
        val = np.cos(np.outer(td, k))                       # f(t_d) = sum_k val[d, k] c_k
        prim = np.empty((nd, n))                            # integral_0^t f = sum_k prim[d, k] c_k
        prim[:, 0] = td
        prim[:, 1:] = np.sin(np.outer(td, k[1:])) / k[1:][None, :]
        _BASIS[key] = (tj, dct, td, val @ dct, prim @ dct)
    return _BASIS[key]


def bb_photons(E, kT):
    return BB_CONST * E * E / np.expm1(E / kT)


class Ring(object):
    """One ring of latitude: cumulative count of the reference element against its arrival angle."""

    def __init__(self, u, f, beta, incl, theta, om_R_over_c, E_obs, kT, tables, n=24, nd=769, beaming=None):
        xt, ct = tables
        self.empty = True
        A = math.cos(incl) * math.cos(theta)
        B = math.sin(incl) * math.sin(theta)
        c_limb = ct[0]
        c_hi, c_lo = A + B, max(A - B, c_limb)
        if c_hi <= c_limb or B < 1e-12:
            return
        self.empty = False
        x_hi = float(np.interp(c_hi, ct, xt))               # two end values per ring (by table here; by a fit later)
        x_lo = float(np.interp(c_lo, ct, xt)) if A - B > c_limb else 0.0
        tj, dct, td, to_dense, to_prim = _basis(n, nd)
        xb, h = 0.5 * (x_hi + x_lo), 0.5 * (x_hi - x_lo)
        x = xb + h * np.cos(tj)
        th, tau = theta_tau(u, np.maximum(x, 1e-12))
        Cx = np.cos(th)
        cphi = np.clip((Cx - A) / B, -1.0, 1.0)
        phi = np.arccos(cphi)
        sphi = np.sqrt(np.maximum(1.0 - cphi * cphi, 1e-300))
        wgeo = h * np.sin(tj) / (B * sphi)
        spsi = np.sin(th)
        salp = np.sqrt(np.maximum(1.0 - x * x, 0.0))
        gam = 1.0 / math.sqrt(1.0 - beta * beta)
        E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
        # On the dense grid the arrival angle is built from cos(psi) and the delay, which are smooth and even in t.
        # (The phase angle itself is not even in t at closest approach, so it must not be put through the cosine series.)
        cphi_d = np.clip((to_dense @ Cx - A) / B, -1.0, 1.0)
        cphi_d[0] = min(1.0, (c_hi - A) / B)
        phi_d = np.arccos(cphi_d)
        tau_d = to_dense @ tau
        self.branches = []
        for s in (-1.0, 1.0):                                # s = -1: before closest approach (approaching side), s = +1: after
            cxi = -s * salp * math.sin(incl) * sphi / spsi
            delta = 1.0 / (gam * (1.0 - beta * cxi))
            mu = x
            Np = bb_photons(E_obs[None, :] / (f * delta[:, None]), kT)
            if beaming is not None:
                Np = Np * beaming(delta * mu)[:, None]
            w = (wgeo * delta**2 / gam * mu)[:, None] * Np   # (n, nE)
            a_d = s * phi_d + om_R_over_c * tau_d            # arrival angle of the reference element, dense in t
            cum_d = to_prim @ w                              # (nd, nE): integral_0^t w
            self.branches.append((a_d, cum_d))
        (am, cm), (ap, cp) = self.branches
        self.tot_minus = cm[-1]                              # counts of the half-turn before closest approach
        self.total = cm[-1] + cp[-1]
        self.a_start = am[-1]                                # arrival angle where the element comes into view (t = pi, minus branch)
        self.a_end = ap[-1]
        # the minus branch runs backwards in t; make arrays increasing in arrival angle
        self.a_all = np.concatenate([am[::-1], ap[1:]])
        self.c_all = np.concatenate([(self.tot_minus[None, :] - cm)[::-1], self.tot_minus[None, :] + cp[1:]], axis=0)

    def cumulative(self, q):
        """Cumulative count at arrival angles q (any real values), shape (len(q), nE)."""
        q = np.asarray(q, dtype=float)
        two_pi = 2.0 * math.pi
        m = np.floor((q - self.a_start) / two_pi)
        qq = q - two_pi * m
        out = np.empty((q.size, self.c_all.shape[1]))
        for e in range(self.c_all.shape[1]):
            out[:, e] = np.interp(qq, self.a_all, self.c_all[:, e], right=self.total[e])
        return out + m[:, None] * self.total[None, :]


def gauss_rings(theta_s, rho, n_per_panel):
    """Latitude nodes and weights across a circular spot, with each ring's half-stretch in longitude.
    Panels are split where the spot starts to cover whole rings (a pole inside the spot); square-root ends are smoothed."""
    g, w = np.polynomial.legendre.leggauss(n_per_panel)
    g = 0.5 * (g + 1.0)
    w = 0.5 * w
    lo, hi = theta_s - rho, theta_s + rho
    panels = []
    if lo < 0.0:                                             # north pole inside the spot
        panels.append((0.0, -lo, False, False))
        panels.append((-lo, min(hi, math.pi), False, True))
    elif hi > math.pi:                                       # south pole inside the spot
        panels.append((lo, 2 * math.pi - hi, True, False))
        panels.append((2 * math.pi - hi, math.pi, False, False))
    else:
        panels.append((lo, hi, True, True))
    th, wt, dl = [], [], []
    for a, b, sa, sb in panels:
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
        t = a + (b - a) * s
        th.append(t)
        wt.append(w * ds * (b - a))
        cd = (math.cos(rho) - np.cos(t) * math.cos(theta_s)) / np.maximum(np.sin(t) * math.sin(theta_s), 1e-300)
        dl.append(np.arccos(np.clip(cd, -1.0, 1.0)))
    return np.concatenate(th), np.concatenate(wt), np.concatenate(dl)


def light_curve(M_sun, R_km, nu_hz, D_pc, incl, theta_s, rho, kT, E_obs, phases, n_rings=12, n=24, beaming=None):
    L = 2.0 * GM_SUN * M_sun / C_CM_S**2
    R = R_km * 1.0e5
    D = D_pc * PC_CM
    u = L / R
    f = math.sqrt(1.0 - u)
    Om = 2.0 * math.pi * nu_hz
    xt = np.concatenate([np.linspace(0.0, 0.02, 400, endpoint=False), np.linspace(0.02, 1.0, 6001)])
    xt[0] = 1e-9
    ct = np.cos(theta_tau(u, xt)[0])
    theta, wth, dlam = gauss_rings(theta_s, rho, n_rings)
    a = np.asarray(phases) * 2.0 * math.pi
    E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
    out = np.zeros((a.size, E_obs.size))
    for th, w, dl in zip(theta, wth, dlam):
        ring = Ring(u, f, Om * R * math.sin(th) / (C_CM_S * f), incl, th, Om * R / C_CM_S, E_obs, kT, (xt, ct), n=n, beaming=beaming)
        if ring.empty:
            continue
        pref = R * R * math.sin(th) * w / (D * D)
        if dl >= math.pi - 1e-12:
            out += pref * ring.total[None, :]
        else:
            out += pref * (ring.cumulative(a + dl) - ring.cumulative(a - dl))
    return out


if __name__ == "__main__":
    from plain_sphere import TESTS
    ref_dir = sys.argv[1]
    for name, (ths, rho, inc, nu) in TESTS.items():
        d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
        y = light_curve(1.4, 12.0, nu, 200.0, math.radians(inc), math.radians(ths), rho, 0.35, 1.0, d[:, 0])[:, 0]
        med = np.median(d[:, 1])
        e = (y - d[:, 1]) / med
        k = int(np.argmax(np.abs(e)))
        big = d[:, 1] > 0.01 * d[:, 1].max()
        print("%s  largest difference over the median %.2e at phase %.4f ; fractional, above 1%% of peak: %.2e ; r.m.s. over the median %.2e"
              % (name, abs(e[k]), d[k, 0], np.max(np.abs(y[big] / d[big, 1] - 1.0)), np.sqrt(np.mean(e**2))), flush=True)
