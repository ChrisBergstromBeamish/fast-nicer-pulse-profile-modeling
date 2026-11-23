# -*- coding: utf-8 -*-
"""Checks of the timed kernel (fast_kernel.py).
1. Against the correctness version (fast_star.py: the pointwise curve, sampled at 64 points inside each of 32 phase bins and
   averaged; it has been compared with the agreed reference curves), on stars of the agreed test families.
2. Against itself at high resolution, to see what its working settings cost in accuracy.
3. Its time per call on this machine, one process."""
import math
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_kernel as fk
import fast_star as fs
from consts import PC_CM

NB, SUB = 32, 64
edges = np.linspace(0.0, 1.0, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
cases = [("SD1e", fs.SD1["SD1e"] + (None,), False, False), ("SD1f", fs.SD1["SD1f"] + (None,), False, False),
         ("OS1b", fs.OS1["OS1b"], True, True), ("OS1e", fs.OS1["OS1e"], True, True), ("OS1f", fs.OS1["OS1f"], True, True),
         ("OS1g", fs.OS1["OS1g"], True, True), ("OS1j", fs.OS1["OS1j"], True, True)]
quick = len(sys.argv) > 1 and sys.argv[1] == "quick"
for name, (ths, rho, inc, nu, bm), ob, ga in ([] if quick else cases):
    ref = fs.light_curve(1.4, 12.0, nu, 200.0, math.radians(inc), math.radians(ths), rho, 0.35, 1.0, mid, n_rings=64,
                         beaming=fs.BEAM[bm], oblate=ob, gamma_area=ga)[:, 0].reshape(NB, SUB).mean(axis=1) / NB
    kw = dict(oblate=ob, gamma_area=ga, beaming=fs.BEAM[bm])
    a = (1.4, 12.0, nu, 200.0 * PC_CM, math.radians(inc), math.radians(ths), rho, 0.35, [1.0], edges)
    hi = fk.spot_counts(*a, n_rings=48, n=32, **kw)[:, 0]
    lo = fk.spot_counts(*a, n_rings=10, n=16, **kw)[:, 0]
    print("%s  kernel at high settings against the correctness version: %.1e ; kernel at working settings against itself at high: %.1e"
          % (name, np.max(np.abs(hi - ref)) / np.mean(ref), np.max(np.abs(lo - hi)) / np.mean(hi)), flush=True)

# the X-PSI quick example's star (314 Hz, 60 deg, spot at 70 deg, 0.75 rad, log10 T = 6.7), 32 bins
K_KEV = 8.617333262e-8
E = np.logspace(-1.0, math.log10(3.0), 48)
args = (1.4, 12.0, 314.0, 1.0e3 * PC_CM, math.radians(60.0), math.radians(70.0), 0.75, 10**6.7 * K_KEV, E, edges)
ref = fk.spot_counts(*args, n_rings=64, n=48)
for nr, n in ((6, 12), (8, 12), (8, 16), (10, 16), (12, 16), (16, 16)):
    S = fk.spot_counts(*args, n_rings=nr, n=n)
    t0 = time.perf_counter()
    for _ in range(30):
        fk.spot_counts(*args, n_rings=nr, n=n)
    dt = (time.perf_counter() - t0) / 30
    print("quick-example star: %2d rings per panel, %2d points: largest bin error over the mean bin %.1e ; time per call %.2f ms"
          % (nr, n, np.max(np.abs(S - ref) / ref.mean(axis=0)), dt * 1e3), flush=True)
