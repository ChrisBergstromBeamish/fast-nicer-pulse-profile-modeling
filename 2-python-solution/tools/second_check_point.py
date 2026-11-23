"""The public code's second shipped check (xpsi/tests/test_likelihood.py): the same data and instrument, no limits on the
background from above, 128 energies, the star p = [1.4, 10, 1, cos 60 deg, 0, 70 deg, 0.75, 6.8]; expected
-47881.27817666349 to a relative tolerance of 1e-5 (0.479). Our programs FAIL this check, by +6.52; the value of that failure
found before (ours converged minus the shipped value) is +6.518, and the RESULT line gives the distance from it.
usage: python -I -B -X utf8 second_check_point.py        (a few seconds; the public example's data from settings.cfg)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
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
RECORDED = 6.518
vals = {}
for nr, n, nE in ((48, 32, 48), (24, 24, 32), (10, 16, 24)):
    v = qe.ours(P2, counts, support, edges_E, n_rings=nr, n=n, nE=nE)
    vals[(nr, n, nE)] = v
    print("ours: %2d rings per panel, %2d points, %2d energies: ln L = %.6f ; minus the public code's expected value: %+.6f (its tolerance 0.479)"
          % (nr, n, nE, v, v - EXP2), flush=True)
star = qe.ours(P2, counts, support, edges_E, 48, 32, 48, return_star=True)
model = qe.T_EXP / qe.NB * star
print("model source counts %.0f against %.0f counts in the data; channels where the model alone exceeds the data: %d of 291"
      % (model.sum(), counts.sum(), int(np.sum(model.sum(axis=1) > counts.sum(axis=1)))))
# sensitivity: change of ln L for a change of the model by one part in 10,000 in every bin
v0 = pl.marginal_loglike(counts, star, qe.T_EXP, support)
v1 = pl.marginal_loglike(counts, star * (1.0 + 1e-4), qe.T_EXP, support)
print("a change of the model counts by 1 part in 10,000 changes ln L by %+.3f at this star" % (v1 - v0))
c0, _, s0 = qe.load_data(settings.data_dir())
st0 = qe.ours(qe.P, c0, s0, edges_E, 48, 32, 48, return_star=True)
w0 = pl.marginal_loglike(c0, st0, qe.T_EXP, s0)
w1 = pl.marginal_loglike(c0, st0 * (1.0 + 1e-4), qe.T_EXP, s0)
print("for comparison, at the quick example's star the same change gives %+.4f" % (w1 - w0))
conv = vals[(48, 32, 48)] - EXP2
print("RESULT converged_minus_shipped=%.6f converged_minus_shipped_vs_recorded=%.3e working_minus_shipped=%.6f sensitivity_second=%.3f sensitivity_quick=%.4f"
      % (conv, abs(conv - RECORDED), vals[(10, 16, 24)] - EXP2, v1 - v0, w1 - w0))
