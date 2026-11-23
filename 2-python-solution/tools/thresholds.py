"""Three numbers the documentation states about L/R = 2GM/(c^2 R): the smallest L/R at which a ray leaving a round star along
its surface sweeps pi (from there a second path of light to a spot at the far pole is possible); the largest L/R the public
example's prior allows (1.6 solar masses, 10 km); and the quick example's star (1.4 solar masses, 12 km). The first is found
by bisection on fast_kernel.bend with 96 Gauss points; the other two use the observed GM of the Sun in consts.py.
usage: python -I -B -X utf8 thresholds.py        (two seconds; no reference material is needed)"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
from consts import C_CM_S, GM_SUN
lo, hi = 0.3, 0.62
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if fk.bend(np.array([mid]), np.array([0.0]), n=96, need_delay=False)[0][0] >= math.pi:
        hi = mid
    else:
        lo = mid
print("L/R at which the ray leaving a round star along its surface sweeps pi: %.4f" % (0.5 * (lo + hi)))
per_sun_per_km = 2.0 * GM_SUN / C_CM_S**2 / 1.0e5
print("largest L/R in the public example's prior (1.6 solar masses, 10 km): %.4f" % (1.6 * per_sun_per_km / 10.0))
print("L/R of the quick example's star (1.4 solar masses, 12 km): %.4f" % (1.4 * per_sun_per_km / 12.0))
