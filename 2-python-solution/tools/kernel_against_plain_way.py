"""The fast kernel against the plain, slow way (plain_sphere.py: a mesh over the spot, the direction search and the lensing
factor by table; a different method) on round stars in many geometries: observers above and below the equator and on the
axis, spots in both halves and over the poles, compact and extended stars. Bin averages over 16 phase bins at 1 keV.
usage: python -I -B -X utf8 kernel_against_plain_way.py        (about a minute and a half; no reference material is needed)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
import math, os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk, plain_sphere as ps
from consts import PC_CM
kT, E1 = 0.35, 1.0
NB, SUB = 16, 16
edges = np.linspace(0, 1, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
rng = np.random.default_rng(11)
cases = [(1.4, 12.0, 400.0, 0.0, 40.0, 0.5), (1.4, 12.0, 400.0, 180.0, 140.0, 0.5), (1.4, 12.0, 400.0, 1e-5, 40.0, 0.5),
         (1.4, 12.0, 300.0, 0.0, 10.0, 0.6), (1.4, 12.0, 300.0, 60.0, 60.0, 0.4), (1.4, 12.0, 300.0, 90.0, 90.0, 1.0),
         (1.48, 11.55, 663.0, 69.0, 156.0, 0.64), (1.62, 10.99, 387.0, 119.0, 36.0, 1.49), (1.26, 11.97, 388.0, 23.0, 153.0, 0.43),
         (1.07, 12.10, 130.0, 74.0, 132.0, 1.08), (1.9, 10.3, 215.0, 88.0, 77.0, 0.34)]
for _ in range(25):
    cases.append((rng.uniform(1.0, 2.0), rng.uniform(10.0, 14.0), rng.uniform(1.0, 700.0), rng.uniform(0.0, 180.0), rng.uniform(3.0, 177.0), rng.uniform(0.05, 1.5)))
worst, nflag, nrefused = 0.0, 0, 0
for (M, R, nu, inc, ths, rho) in cases:
    i, t = math.radians(inc), math.radians(ths)
    u = 2.953250 * M / R
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c = fk.spot_counts(M, R, nu, 200.0 * PC_CM, i, t, rho, kT, [E1], edges, n_rings=24, n=24, oblate=False, gamma_area=False)[:, 0]
    except (fk.TwoPaths, ValueError) as ex:
        nrefused += 1
        print("refused (L/R = %.3f): %s" % (u, str(ex)[:90]))
        continue
    a = ps.light_curve(M, R, nu, 200.0, i, t, rho, kT, E1, mid, n_zeta=48, n_eta=192, n_phi=2048).reshape(NB, SUB).mean(axis=1) / NB
    sc = a.max() if a.max() > 0 else 1.0
    e = np.max(np.abs(c - a)) / sc
    ok = np.all(np.isfinite(c)) and (e < 2e-4 or a.max() < 1e-12)
    worst = max(worst, e) if np.isfinite(e) else float("inf")
    if not ok:
        nflag += 1
    print("%s observer at %5.1f deg, spot at %5.1f deg, radius %.2f, L/R %.3f, %3.0f Hz: plain mean bin %.3e ; kernel minus plain %.1e of the largest bin"
          % ("  " if ok else "**", inc, ths, rho, u, nu, a.mean(), e))
print("%d geometries: %d flagged (difference above 2e-4), %d refused, largest difference %.1e" % (len(cases), nflag, nrefused, worst))
print("RESULT geometries=%d flagged=%d refused=%d largest=%.3e" % (len(cases), nflag, nrefused, worst))
