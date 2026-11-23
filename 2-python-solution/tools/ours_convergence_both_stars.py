# -*- coding: utf-8 -*-
"""Is "ours, converged" (48 rings per panel, 32 points, 48 energies) converged at both of the public code's check stars,
the quick example's and the second check's? Ours at finer settings beside it: the value of ln L and the bins (total, the
largest single bin, the largest per phase bin summed over channels, the largest per channel summed over the turn), each
against the 48-32-48 setting. Only our own programs are run. Written 23 Nov 2026 for the comparison with the public code.
usage: python -I -B -X utf8 2-python-solution/tools/ours_convergence_both_stars.py        (the public example's data folder
       comes from settings.cfg)"""
import math
import os
import sys
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_v] = "1"          # one thread, set before NumPy loads
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk
import public_likelihood as pl
import quick_example as qe
import settings
from consts import BB_CONST

counts, bg, support = qe.load_data(settings.data_dir())
nosup = np.stack([np.zeros(counts.shape[0]), -np.ones(counts.shape[0])], axis=1)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph_edges = np.linspace(0.0, 1.0, qe.NB + 1)


def star_counts(p, nr, n, nE):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W
    mass, radius, distance, cos_i, shift, colat, rho, logT = p
    M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
    S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                       qe.X_KB_OVER_KEV * 10.0**logT, E, ph_edges, phase0=0.5 - shift, n_rings=nr, n=n)
    return Wn @ S.T


for name, p, sup in (("quick", [1.4, 12.0, 1.0, math.cos(60 * math.pi / 180), 0.0, 70 * math.pi / 180, 0.75, 6.7], support),
                     ("second", [1.4, 10.0, 1.0, math.cos(60 * math.pi / 180), 0.0, 70 * math.pi / 180, 0.75, 6.8], nosup)):
    ref = star_counts(p, 48, 32, 48)
    v_ref = pl.marginal_loglike_fast(counts, ref, qe.T_EXP, sup, fact)
    print("%s check: ours at 48 rings per panel, 32 points, 48 energies: ln L = %.6f" % (name, v_ref))
    for nr, n, nE in ((10, 16, 24), (48, 32, 96), (48, 64, 48), (96, 32, 48), (96, 64, 96), (144, 96, 144)):
        s = star_counts(p, nr, n, nE)
        v = pl.marginal_loglike_fast(counts, s, qe.T_EXP, sup, fact)
        rel = s / ref - 1.0
        print("   %3d rings, %2d points, %3d energies: ln L = %.6f ; minus the one above %+.6f ; bins: total %+.1e ; largest %.1e ; per phase bin largest %.1e ; per channel largest %.1e"
              % (nr, n, nE, v, v - v_ref, s.sum() / ref.sum() - 1.0, np.max(np.abs(rel)),
                 np.max(np.abs(s.sum(axis=0) / ref.sum(axis=0) - 1.0)), np.max(np.abs(s.sum(axis=1) / ref.sum(axis=1) - 1.0))), flush=True)
