"""The even grid that locates arrival times (fast_kernel.ND, 129 points at the working setting) varied on its own, from 65 to
2049 points, on the quick example at the working setting (10 rings per panel, 16 points, 24 energies): how much the bins and the
likelihood move against the finest grid. The point where the arrival angle equals a bin edge is found on this grid and not
refined; the fast_kernel.py docstring says what the measured effect is.
usage: python -I -B -X utf8 arrival_grid_check.py        (a few seconds; the public example's data from settings.cfg)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk, quick_example as qe, public_likelihood as pl, settings
from consts import BB_CONST
counts, bg, support = qe.load_data(settings.data_dir())
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], 24, edges_E)
M, D, kT, inc = qe.star_args(qe.P)
ph = np.linspace(0.0, 1.0, qe.NB + 1)
grids = (65, 129, 257, 513, 1025, 2049)
out = {}
for nd in grids:
    fk.ND = nd
    fk._CH.clear()                      # the Chebyshev tools hold the grid; made again for each size
    S = fk.spot_counts(M, qe.P["radius"], qe.P["freq"], D, inc, qe.P["colat"], qe.P["ang_radius"], kT, E, ph,
                       phase0=0.5 - qe.P["phase_shift"], n_rings=10, n=16) * (qe.X_BB / BB_CONST)
    v = pl.marginal_loglike_fast(counts, (qe.AREA * qe.NB) * (W @ S.T), qe.T_EXP, support, fact)
    out[nd] = (S, v)
fk.ND = 129
fk._CH.clear()
S_ref, v_ref = out[2049]
res = {}
for nd in grids:
    S, v = out[nd]
    eb = float(np.max(np.abs(S - S_ref)) / S_ref.mean())
    print("grid of %4d points: ln L = %.9f ; minus the 2049-point grid %+.2e ; bins differ by %.1e of the mean bin" % (nd, v, v - v_ref, eb))
    res[nd] = (abs(v - v_ref), eb)
print("RESULT lnL_129=%.3e bins_129=%.3e lnL_65=%.3e bins_65=%.3e" % (res[129][0], res[129][1], res[65][0], res[65][1]))
