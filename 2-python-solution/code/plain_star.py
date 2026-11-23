# -*- coding: utf-8 -*-
"""The usual way, kept plain and slow, for a round or a flattened star: a mesh of small patches over the spot, a fine grid of
emission phases, for each patch its own table of light's path (direction against angle, the lensing factor by differences of
that table, the delay), the flux of each patch resampled at the arrival phases, and the patches added.

It shares with the fast kernel only the path integrals (fast_kernel.bend with its own number of points) and the inputs.
The method is different at every other step: no count at emission, no rings, no exact stretches, no function K.
It exists to check the fast kernel bin by bin (Codex's round-3 test: an independent route in emission phase with the
lensing factor written out).

Photon flux from a patch (NICER_SPEED_PLAN_v0.2.md section 3):  dF = (dS'/D^2) J delta^3 mu N'(E/(f delta)),  J = |d cos(alpha)/d cos(psi)|,
dS' = gamma R^2 sqrt(1 + h^2) dOmega for the flattened family (the spot is a cap in the star's own angles), dOmega the patch's
solid angle seen from the centre."""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_kernel as fk
from consts import BB_CONST, C_CM_S, GM_SUN


def light_curve(M_sun, Req_km, nu_hz, D_cm, incl, theta_s, rho, kT, E_obs, phases, n_zeta=48, n_eta=192, n_phi=2048, n_tab=400,
                oblate=True, gamma_area=True, n_gl=48, chunk=128):
    """Flux at the given phases (turns; the spot's centre is closest at phase 0) and energies: array (phases, energies)."""
    GM = GM_SUN * M_sun
    L = 2.0 * GM / C_CM_S**2
    Req = Req_km * 1.0e5
    Om = 2.0 * math.pi * nu_hz
    o2 = (Om**2 * Req**3 / GM) * (-0.788 + 1.030 * GM / (Req * C_CM_S**2)) if oblate else 0.0
    E_obs = np.atleast_1d(np.asarray(E_obs, dtype=float))
    # mesh over the spot in its own polar angles; each patch's solid angle
    ze = np.linspace(0.0, rho, n_zeta + 1)
    zc = np.arccos(0.5 * (np.cos(ze[:-1]) + np.cos(ze[1:])))
    dOm = (np.cos(ze[:-1]) - np.cos(ze[1:])) * (2.0 * math.pi / n_eta)
    eta = (np.arange(n_eta) + 0.5) * 2.0 * math.pi / n_eta
    Z, ETA = np.meshgrid(zc, eta, indexing="ij")
    ns = np.array([math.sin(theta_s), 0.0, math.cos(theta_s)])
    et = np.array([math.cos(theta_s), 0.0, -math.sin(theta_s)])
    px = np.cos(Z) * ns[0] + np.sin(Z) * np.cos(ETA) * et[0]
    py = np.sin(Z) * np.sin(ETA)
    pz = np.cos(Z) * ns[2] + np.sin(Z) * np.cos(ETA) * et[2]
    theta = np.arccos(np.clip(pz, -1.0, 1.0)).ravel()
    lam = np.arctan2(py, px).ravel()
    dOm = np.repeat(dOm[:, None], n_eta, axis=1).ravel()
    ct_, st_ = np.cos(theta), np.sin(theta)
    R = Req * (1.0 + o2 * ct_ * ct_)
    u = L / R
    f = np.sqrt(1.0 - u)
    h = -2.0 * Req * o2 * ct_ * st_ / (R * f)
    rt = np.sqrt(1.0 + h * h)
    cg, sg = 1.0 / rt, h / rt
    beta = Om * R * st_ / (C_CM_S * f)
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    area = R * R * rt * dOm * (gam if gamma_area else 1.0)
    trad = (Req - R + L * np.log((Req - L) / (R - L))) / C_CM_S
    ci, si = math.cos(incl), math.sin(incl)
    phi_e = (np.arange(n_phi) + 0.5) * 2.0 * math.pi / n_phi
    two_pi = 2.0 * math.pi
    ph_out = np.asarray(phases, dtype=float) * two_pi
    out = np.zeros((ph_out.size, E_obs.size))
    # each patch's own table: x from its floor (a little inward of the tilt) to 1, denser near the floor
    s_tab = np.linspace(0.0, 1.0, n_tab)
    for c0 in range(0, theta.size, chunk):
        sl = slice(c0, min(c0 + chunk, theta.size))
        x_floor = -(np.abs(sg[sl]) + 0.02)
        xt = x_floor[:, None] + (1.0 - x_floor)[:, None] * (s_tab ** 1.5)[None, :]         # (cells, n_tab), increasing
        xt[:, -1] = 1.0 - 1e-13
        th, ta = fk.bend(u[sl][:, None], xt, n=n_gl)
        if not np.all(th < math.pi):
            raise fk.TwoPaths("the plain way also adds one path only")
        ctab = np.cos(th)                                                                    # increasing with x
        for m, k in enumerate(range(sl.start, sl.stop)):
            pc = lam[k] + phi_e
            cpsi = ci * ct_[k] + si * st_[k] * np.cos(pc)
            reach = cpsi > ctab[m, 0]
            x = np.interp(cpsi, ctab[m], xt[m])
            J = np.interp(cpsi, ctab[m], np.gradient(xt[m], ctab[m]))
            tau = np.interp(cpsi, ctab[m], ta[m])
            spsi = np.sqrt(np.maximum(1.0 - cpsi * cpsi, 1e-300))
            salpha = np.sqrt(np.maximum(1.0 - x * x, 0.0))
            mu = x * cg[k] + salpha * sg[k] * (ci - ct_[k] * cpsi) / (st_[k] * spsi)
            cxi = -salpha * si * np.sin(pc) / spsi
            delta = 1.0 / (gam[k] * (1.0 - beta[k] * cxi))
            vis = reach & (mu > 0.0)
            geom = np.where(vis, (area[k] / (D_cm * D_cm)) * np.abs(J) * delta**3 * mu, 0.0)
            a = phi_e + Om * (tau * R[k] / C_CM_S + trad[k])
            ae = np.concatenate([a - two_pi, a, a + two_pi])
            Ep = E_obs[None, :] / (f[k] * delta[:, None])
            F = geom[:, None] * (BB_CONST * Ep * Ep / np.expm1(Ep / kT))
            for e in range(E_obs.size):
                out[:, e] += np.interp(ph_out, ae, np.concatenate([F[:, e], F[:, e], F[:, e]]))
    return out
