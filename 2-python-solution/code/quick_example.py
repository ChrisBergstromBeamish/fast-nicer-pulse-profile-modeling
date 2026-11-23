# -*- coding: utf-8 -*-
"""The public code's shipped quick example (X-PSI, examples/examples_fast/Modules/main.py), worked with our own programs.

The example: one star (the vector p), one data set (291 channels, 32 phase bins, 1000 s), and the log-likelihood the public
code expects for it: -3.2552616905e+04, which its own check accepts to a relative tolerance of 1e-5 (0.3255).
X-PSI is not imported and not run. Its data files are read as numbers; its constants are used as the example's inputs.

Two ways are worked here:
  ours      bin counts straight from fast_kernel.spot_counts (exact phase-bin integrals), energy integrals by a Chebyshev
            fit in ln(E), then the public code's likelihood definition (public_likelihood.py).
  emulated  our converged curve sampled at the public code's 64 phases and 64 energies, then the public code's own
            interpolation steps (akima-periodic in energy and in phase). It shows what those two steps cost in accuracy.
"""
import math
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_kernel as fk
import fast_star as fs
import public_likelihood as pl
from consts import BB_CONST, C_CM_S, GM_SUN, PC_CM

# ---- the public code's constants (xpsi/global_imports.py and hot_BB.pyx); inputs of the example, not ours
X_GM_M = 0.5 * 2.95325024e3                    # G M_sun / c^2 in metres
X_KPC_M = 3.08567758e19
X_KB_OVER_KEV = 1.38064852e-23 / 1.60217662e-16
X_BB = 5.040366110812353e22 * 1.0e-7 / 1.60217662e-16     # photons / cm^2 / s / keV / sr per keV^2
EXPECTED = -3.2552616905e+04
TOL = 1.0e-5 * abs(EXPECTED)

P = dict(mass=1.4, radius=12.0, distance=1.0, cos_i=math.cos(60.0 * math.pi / 180.0), phase_shift=0.0,
         colat=70.0 * math.pi / 180.0, ang_radius=0.75, log10T=6.7, freq=314.0)
T_EXP = 1000.0
AREA = 1800.0
NB = 32


def star_args(P):
    M_eff = P["mass"] * X_GM_M * 100.0 * C_CM_S**2 / GM_SUN      # so that G M / c^2 is the public code's value
    D_cm = P["distance"] * X_KPC_M * 100.0
    kT = X_KB_OVER_KEV * 10.0 ** P["log10T"]
    return M_eff, D_cm, kT, math.acos(P["cos_i"])


def load_data(data_dir):
    counts = np.loadtxt(os.path.join(data_dir, "xpsi_good_realisation.dat"), dtype=float)
    bg = np.loadtxt(os.path.join(data_dir, "xpsi_synthetic_background.dat"))[:, 1]
    support = np.stack([0.9 * bg / T_EXP, 1.1 * bg / T_EXP], axis=1)
    return counts, bg, support


def cheb_energy_matrix(lo, hi, N, edges):
    """Nodes E_k (Chebyshev points in ln E on [lo, hi]) and the matrix W with
    INT_{edges[q]}^{edges[q+1]} F(E) dE  ~  sum_k W[q, k] F(E_k)   (a polynomial of degree N-1 in ln E through E F)."""
    a, b = math.log(lo), math.log(hi)
    k = np.arange(N)
    thj = math.pi * (k + 0.5) / N
    x = 0.5 * (a + b) + 0.5 * (b - a) * np.cos(thj)
    V = np.cos(np.outer(k, thj)) * (2.0 / N)
    V[0] *= 0.5
    xi = np.clip((2.0 * np.log(edges) - (a + b)) / (b - a), -1.0, 1.0)
    T = np.cos(np.outer(np.arccos(xi), np.arange(N + 1)))                 # T_m(xi), m = 0 .. N
    Pm = np.zeros((xi.size, N))
    Pm[:, 0] = T[:, 1]
    if N > 1:
        Pm[:, 1] = 0.25 * T[:, 2]
    for m in range(2, N):
        Pm[:, m] = 0.5 * (T[:, m + 1] / (m + 1) - T[:, m - 1] / (m - 1))
    W = 0.5 * (b - a) * ((Pm[1:] - Pm[:-1]) @ V)                          # for g(x) = E F(E) at the nodes
    E = np.exp(x)
    return E, W * E[None, :]


def ours(P, counts, support, edges_E, n_rings=10, n=16, nE=24, return_star=False):
    M_eff, D_cm, kT, incl = star_args(P)
    E, W = cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    ph_edges = np.linspace(0.0, 1.0, NB + 1)
    S = fk.spot_counts(M_eff, P["radius"], P["freq"], D_cm, incl, P["colat"], P["ang_radius"], kT, E, ph_edges,
                       phase0=0.5 - P["phase_shift"], n_rings=n_rings, n=n) * (X_BB / BB_CONST)
    star = (AREA * NB) * (W @ S.T)                                        # (channel, bin): mean count rate in the bin
    if return_star:
        return star
    return pl.marginal_loglike(counts, star, T_EXP, support)


def emulated(P, counts, support, edges_E, n_rings=64, n=24, nd=769, return_star=False):
    M_eff, D_cm, kT, incl = star_args(P)
    energies = pl.example_energies(edges_E, 64)
    phases = np.linspace(0.0, 1.0, 64)
    F = pointwise(P, energies, phases, n_rings=n_rings, n=n, nd=nd)       # (phase, energy)
    binned = pl.energy_integrate(F.T, energies, edges_E)                  # (interval, phase)
    rate = AREA * binned
    star = pl.phase_integrate(rate, phases, np.linspace(0.0, 1.0, NB + 1), P["phase_shift"])
    if return_star:
        return star
    return pl.marginal_loglike(counts, star, T_EXP, support)


def pointwise(P, energies, phases, n_rings=64, n=24, nd=769):
    """Our curve at given phases (cycles, the public code's phase: the spot's centre is opposite the observer at phase 0)."""
    M_eff, D_cm, kT, incl = star_args(P)
    star = fs.Star(M_eff, P["radius"], P["freq"], oblate=True)
    theta, wth, dlam = fk.spot_rings(P["colat"], P["ang_radius"], n_rings,
                                     breaks=fk.visibility_breaks(star.L, star.Req, star.o2, incl))
    a = (np.asarray(phases) - 0.5) * 2.0 * math.pi
    out = np.zeros((a.size, energies.size))
    for th, w, dl in zip(theta, wth, dlam):
        ring = fs.Ring(star, incl, th, energies, kT, n=n, nd=nd)
        if ring.empty:
            continue
        pref = ring.area * w / (D_cm * D_cm)
        if dl >= math.pi - 1e-12:
            out += pref * ring.total[None, :]
        else:
            out += pref * (ring.cumulative(a + dl) - ring.cumulative(a - dl))
    return out * (X_BB / BB_CONST)


if __name__ == "__main__":
    # usage: python -I -B -X utf8 quick_example.py [emulated]   (the public example's data folder comes from settings.cfg)
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
    import settings
    data_dir = settings.data_dir()
    counts, bg, support = load_data(data_dir)
    edges_E = pl.example_instrument_edges()
    print("data: %d channels x %d bins, %d counts ; instrument intervals %.4f to %.4f keV ; expected %.6f, tolerance %.4f"
          % (counts.shape[0], counts.shape[1], counts.sum(), edges_E[0], edges_E[-1], EXPECTED, TOL), flush=True)
    for nr, n, nE in ((48, 32, 48), (24, 24, 32), (16, 16, 24), (12, 16, 24), (10, 16, 24), (8, 16, 24), (8, 12, 16), (6, 12, 16)):
        v = ours(P, counts, support, edges_E, n_rings=nr, n=n, nE=nE)
        print("ours: %2d rings per panel, %2d points, %2d energies: ln L = %.6f ; minus the public code's expected value: %+.6f"
              % (nr, n, nE, v, v - EXPECTED), flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == "emulated":
        for nr, n, nd in ((64, 24, 769), (64, 32, 3073)):
            t0 = time.perf_counter()
            v = emulated(P, counts, support, edges_E, n_rings=nr, n=n, nd=nd)
            print("emulated steps of the public code on our curve (%d rings per panel, %d points, grid %d): ln L = %.6f ; minus expected: %+.6f  (%.0f s)"
                  % (nr, n, nd, v, v - EXPECTED, time.perf_counter() - t0), flush=True)
