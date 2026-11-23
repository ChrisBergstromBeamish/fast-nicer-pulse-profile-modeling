"""The public code's two interpolation steps (its 128 energies, 64 phases) applied to our curve at the public code's second check
point (radius 10 km, log10 T 6.8, no upper limit on the background).
usage: python -I -B -X utf8 second_check_point_emulated.py        (several minutes; the public example's data folder from settings.cfg)"""
import math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quick_example as qe
import public_likelihood as pl
import settings
counts, bg, _ = qe.load_data(settings.data_dir())
support = np.stack([np.zeros(291), -np.ones(291)], axis=1)
edges_E = pl.example_instrument_edges()
P2 = dict(qe.P, radius=10.0, log10T=6.8)
EXP2 = -47881.27817666349
phases = np.linspace(0.0, 1.0, 64)
for nE in (128,):
    energies = pl.example_energies(edges_E, nE)
    for nr, n, nd in ((64, 24, 769), (64, 32, 3073), (96, 32, 6145)):
        F = qe.pointwise(P2, energies, phases, n_rings=nr, n=n, nd=nd)
        rate = qe.AREA * pl.energy_integrate(F.T, energies, edges_E)
        star = pl.phase_integrate(rate, phases, np.linspace(0.0, 1.0, 33), 0.0)
        v = pl.marginal_loglike(counts, star, qe.T_EXP, support)
        print("the public code's two interpolation steps (%d energies, 64 phases) on our curve (%d rings per panel, %d points, grid %d): ln L = %.6f ; minus expected: %+.6f"
              % (nE, nr, n, nd, v, v - EXP2), flush=True)
