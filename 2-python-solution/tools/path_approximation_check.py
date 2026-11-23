"""The approximation of the path integrals (fast_kernel._path_fit: one two-variable Chebyshev approximation of the swept angle
and the delay for every star in its range) against the integrals themselves (fast_kernel.bend with 96 Gauss points), on 6000
random directions: u = L/R over the approximation's range, x = cos(alpha) from its lowest value to 1.
usage: python -I -B -X utf8 path_approximation_check.py        (a few seconds; no reference material is needed)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
rng = np.random.default_rng(5)
N = 6000
u = rng.uniform(fk.FIT_U[0], fk.FIT_U[1], N)
x = rng.uniform(fk.FIT_XMIN, 1.0, N)
th_int, ta_int = fk.bend(u, x, n=96)
cth, cta, ZM, kz = fk._ring_series(u)
th_fit, ta_fit = fk._series_eval(cth, cta, ZM, kz, x[:, None])
e_th = float(np.max(np.abs(th_fit[:, 0] - th_int)))
e_ta = float(np.max(np.abs(ta_fit[:, 0] - ta_int)))
print("%d random directions, u in [%.2f, %.2f], x in [%.2f, 1]: the approximation minus the integrals (96 points): swept angle largest %.1e rad ; delay largest %.1e R/c"
      % (N, fk.FIT_U[0], fk.FIT_U[1], fk.FIT_XMIN, e_th, e_ta))
# the integrals against themselves with more points, on the same directions
th_more, ta_more = fk.bend(u, x, n=160)
print("the integrals with 96 points against 160 points: swept angle %.1e rad ; delay %.1e R/c" % (np.max(np.abs(th_more - th_int)), np.max(np.abs(ta_more - ta_int))))
print("RESULT directions=%d swept_angle=%.3e delay=%.3e integrals_96_vs_160_swept=%.3e integrals_96_vs_160_delay=%.3e"
      % (N, e_th, e_ta, np.max(np.abs(th_more - th_int)), np.max(np.abs(ta_more - ta_int))))
