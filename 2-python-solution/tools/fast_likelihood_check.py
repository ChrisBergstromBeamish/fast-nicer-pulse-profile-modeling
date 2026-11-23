"""The all-channels-at-once likelihood against the reference form, on both of the public code's check points, and timed.
usage: python -I -B -X utf8 fast_likelihood_check.py        (the public example's data folder comes from settings.cfg)"""
import math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quick_example as qe
import public_likelihood as pl
import settings
counts, bg, support = qe.load_data(settings.data_dir())
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
nosup = np.stack([np.zeros(291), -np.ones(291)], axis=1)
wide = np.stack([0.5 * support[:, 0], 3.0 * support[:, 1]], axis=1)
P2 = dict(qe.P, radius=10.0, log10T=6.8)
P3 = dict(qe.P, radius=12.6, log10T=6.66, distance=1.3)
for label, P, sup in (("quick example's star, background limited to 0.9 to 1.1 of the given one", qe.P, support),
                      ("second check point, no upper limit on the background", P2, nosup),
                      ("quick example's star, no upper limit on the background", qe.P, nosup),
                      ("quick example's star, background limited to 0.45 to 3.3 of the given one", qe.P, wide),
                      ("a fainter star (12.6 km, 10^6.66 K, 1.3 kpc), background limited to 0.9 to 1.1", P3, support)):
    star = qe.ours(P, counts, sup, edges_E, 10, 16, 24, return_star=True)
    ref = pl.marginal_loglike(counts, star, qe.T_EXP, sup)
    det = {}
    fast = pl.marginal_loglike_fast(counts, star, qe.T_EXP, sup, fact, details=det)
    best = 1e9
    for _ in range(5):
        t0 = time.perf_counter()
        for _ in range(100):
            pl.marginal_loglike_fast(counts, star, qe.T_EXP, sup, fact)
        best = min(best, (time.perf_counter() - t0) / 100)
    f = det["fall"]
    cls = np.bincount(np.digitize(f[f >= 0.0], [0.5, 2.0, 8.0, 18.0, 45.0]), minlength=6)
    print("%s:" % label)
    print("   reference %.9f ; fast %.9f ; difference %.2e ; fast form %.3f ms per call" % (ref, fast, fast - ref, best * 1e3))
    print("   channels by the power series: %d ; the others by number of points (10,14,20,28,40,64): %s"
          % (int(det["series"].sum()), cls), flush=True)
