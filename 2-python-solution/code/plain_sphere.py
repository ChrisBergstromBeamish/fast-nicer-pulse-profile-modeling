# -*- coding: utf-8 -*-
"""The usual way, kept plain and slow, for a round star (the test article's SD1 problems): a mesh over the spot, a fine grid of
emission phases, the direction search by table, the lensing factor by table differences, and resampling at the arrival phases.
It exists to check the formulas and conventions against the agreed reference curves, and as an independent check of the fast code.
Photon flux arriving per unit time, energy and detector area from a surface element (derivation: NICER_SPEED_PLAN_v0.2.md section 3):
    dF = (dS'/D^2) J delta^3 mu N'(E/(f delta); delta mu),   J = d cos(alpha)/d cos(psi),
with dS' the element's area in the rotating star's own frame (the test problems give the spot in that frame)."""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bend import theta_tau
from consts import BB_CONST, C_CM_S, GM_SUN, PC_CM


def bb_photons(E, kT):
    """Black-body photon intensity in the surface's own frame, photons/(cm^2 s sr keV); isotropic."""
    return BB_CONST * E * E / np.expm1(E / kT)


def light_curve(M_sun, R_km, nu_hz, D_pc, incl, theta_s, rho, kT, E_obs, phases, n_zeta=24, n_eta=96, n_phi=2048, beaming=None):
    L = 2.0 * GM_SUN * M_sun / C_CM_S**2
    R = R_km * 1.0e5
    D = D_pc * PC_CM
    u = L / R
    f = math.sqrt(1.0 - u)
    Om = 2.0 * math.pi * nu_hz
    # tables of the forward function on a fine grid of cos(alpha)
    xt = np.concatenate([np.linspace(0.0, 0.02, 400, endpoint=False), np.linspace(0.02, 1.0, 6001)])
    xt[0] = 1e-9
    th, ta = theta_tau(u, xt)
    ct = np.cos(th)                      # increasing with x
    jt = np.gradient(xt, ct)             # lensing factor
    c_limb = ct[0]
    # mesh over the spot in its own polar angles
    ze = np.linspace(0.0, rho, n_zeta + 1)
    zc = np.arccos(0.5 * (np.cos(ze[:-1]) + np.cos(ze[1:])))
    dA = (np.cos(ze[:-1]) - np.cos(ze[1:])) * (2.0 * math.pi / n_eta) * R * R
    eta = (np.arange(n_eta) + 0.5) * 2.0 * math.pi / n_eta
    Z, ETA = np.meshgrid(zc, eta, indexing="ij")
    AREA = np.repeat(dA[:, None], n_eta, axis=1)
    ns = np.array([math.sin(theta_s), 0.0, math.cos(theta_s)])
    et = np.array([math.cos(theta_s), 0.0, -math.sin(theta_s)])
    el = np.array([0.0, 1.0, 0.0])
    px = np.cos(Z) * ns[0] + np.sin(Z) * (np.cos(ETA) * et[0] + np.sin(ETA) * el[0])
    py = np.cos(Z) * ns[1] + np.sin(Z) * (np.cos(ETA) * et[1] + np.sin(ETA) * el[1])
    pz = np.cos(Z) * ns[2] + np.sin(Z) * (np.cos(ETA) * et[2] + np.sin(ETA) * el[2])
    theta = np.arccos(np.clip(pz, -1.0, 1.0)).ravel()
    lam = np.arctan2(py, px).ravel()
    area = AREA.ravel()
    phi_e = (np.arange(n_phi) + 0.5) * 2.0 * math.pi / n_phi
    out = np.zeros(len(phases))
    two_pi = 2.0 * math.pi
    ph_out = np.asarray(phases) * two_pi
    for k in range(theta.size):
        st, ctk = math.sin(theta[k]), math.cos(theta[k])
        pc = lam[k] + phi_e                                   # the cell's angle from the observer's meridian
        cpsi = math.cos(incl) * ctk + math.sin(incl) * st * np.cos(pc)
        vis = cpsi > c_limb
        x = np.interp(cpsi, ct, xt)
        J = np.interp(cpsi, ct, jt)
        tau = np.interp(cpsi, ct, ta)
        spsi = np.sqrt(np.maximum(1.0 - cpsi * cpsi, 1e-300))
        salpha = np.sqrt(np.maximum(1.0 - x * x, 0.0))
        beta = Om * R * st / (C_CM_S * f)
        gam = 1.0 / math.sqrt(1.0 - beta * beta)
        cxi = -salpha * math.sin(incl) * np.sin(pc) / spsi
        delta = 1.0 / (gam * (1.0 - beta * cxi))
        mu = x
        Eprime = E_obs / (f * delta)
        Np = bb_photons(Eprime, kT)
        if beaming is not None:
            Np = Np * beaming(delta * mu)
        # the test problems give the spot in the rotating star's own angles, so area[k] is already the surface's own area dS'
        F = np.where(vis, (area[k] / (D * D)) * J * delta**3 * mu * Np, 0.0)
        a = phi_e + Om * tau * R / C_CM_S                      # arrival phase (radians), monotonic in phi_e
        ae = np.concatenate([a - two_pi, a, a + two_pi])
        Fe = np.concatenate([F, F, F])
        out += np.interp(ph_out, ae, Fe)
    return out


TESTS = {  # spot colatitude (deg), spot angular radius (rad), observer colatitude (deg), spin (Hz)
    "SD1a": (90.0, 0.01, 90.0, 1.0), "SD1b": (90.0, 1.0, 90.0, 1.0), "SD1c": (90.0, 0.01, 90.0, 200.0),
    "SD1d": (90.0, 1.0, 90.0, 200.0), "SD1e": (60.0, 1.0, 30.0, 400.0), "SD1f": (20.0, 1.0, 80.0, 400.0)}

if __name__ == "__main__":
    ref_dir = sys.argv[1]
    for name, (ths, rho, inc, nu) in TESTS.items():
        d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
        small = rho < 0.1
        y = light_curve(1.4, 12.0, nu, 200.0, math.radians(inc), math.radians(ths), rho, 0.35, 1.0, d[:, 0],
                        n_zeta=3 if small else 40, n_eta=12 if small else 160, n_phi=4096)
        med = np.median(d[:, 1])
        big = d[:, 1] > 0.01 * d[:, 1].max()
        print("%s  largest difference over the median %.2e ; largest fractional difference where the flux is above 1%% of its peak %.2e"
              % (name, np.max(np.abs(y - d[:, 1])) / med, np.max(np.abs(y[big] / d[big, 1] - 1.0))), flush=True)
