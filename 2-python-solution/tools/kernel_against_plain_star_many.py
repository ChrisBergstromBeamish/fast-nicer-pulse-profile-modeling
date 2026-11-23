"""The fast kernel against the plain way for flattened stars in many geometries, including fast spins (strong flattening),
observers below the equator and on the axis, spots over the poles. 16 phase bins, two energies. The plain way: 64 x 256
patches, 32 samples inside each bin."""
import math, os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk, plain_star as pst
from consts import PC_CM
kT = 0.35
E = np.array([0.3, 1.0])
NB, SUB = 16, 32
edges = np.linspace(0, 1, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
rng = np.random.default_rng(31)
cases = [(1.4, 12.0, 314.0, 60.0, 70.0, 0.75), (1.4, 10.0, 314.0, 60.0, 70.0, 0.75), (1.4, 12.0, 600.0, 30.0, 60.0, 1.0), (1.4, 12.0, 600.0, 80.0, 20.0, 1.0),
         (1.4, 12.0, 600.0, 0.0, 50.0, 0.6), (1.4, 12.0, 600.0, 43.1, 139.2, 0.3), (1.6, 13.0, 716.0, 55.8, 129.7, 0.4), (1.2, 12.5, 716.0, 114.4, 62.1, 0.5)]
for _ in range(16):
    cases.append((rng.uniform(1.0, 2.0), rng.uniform(10.0, 13.5), float(rng.choice([200.0, 314.0, 500.0, 600.0])), rng.uniform(0.0, 180.0), rng.uniform(3.0, 177.0), rng.uniform(0.05, 1.5)))
worst, nflag, nref = 0.0, 0, 0
for (M, R, nu, inc, ths, rho) in cases:
    i, t = math.radians(inc), math.radians(ths)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c = fk.spot_counts(M, R, nu, 200.0 * PC_CM, i, t, rho, kT, E, edges, n_rings=24, n=24)
            c10 = fk.spot_counts(M, R, nu, 200.0 * PC_CM, i, t, rho, kT, E, edges, n_rings=10, n=16)
            a = pst.light_curve(M, R, nu, 200.0 * PC_CM, i, t, rho, kT, E, mid, n_zeta=64, n_eta=256, n_phi=2048).reshape(NB, SUB, 2).mean(axis=1) / NB
    except (fk.TwoPaths, ValueError, ArithmeticError) as ex:
        nref += 1
        print("refused: %.2f Suns, %.2f km, %3.0f Hz, observer %5.1f, spot %5.1f, radius %.2f: %s" % (M, R, nu, inc, ths, rho, str(ex)[:80]), flush=True)
        continue
    sc = a.max() if a.max() > 0 else 1.0
    e, e10 = np.max(np.abs(c - a)) / sc, np.max(np.abs(c10 - a)) / sc
    ok = np.all(np.isfinite(c)) and e < 2e-4
    worst = max(worst, e) if np.isfinite(e) else float("inf")
    nflag += 0 if ok else 1
    print("%s %.2f Suns, %5.2f km, %3.0f Hz, observer at %5.1f deg, spot at %5.1f deg, radius %.2f: plain mean bin %.2e ; kernel (24 rings, 24 points) minus plain %.1e ; kernel (10, 16) minus plain %.1e"
          % ("  " if ok else "**", M, R, nu, inc, ths, rho, a.mean(), e, e10), flush=True)
print("%d geometries: %d flagged (above 2e-4 of the largest bin), %d refused, largest difference %.1e" % (len(cases), nflag, nref, worst))
