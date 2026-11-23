"""The fast kernel against the plain way for flattened stars (plain_star.py), bin by bin and energy by energy, on the public
quick example's star and on its second check star (an independent route in emission phase with the lensing factor written
out). Also the plain way itself against three of the agreed reference curves. The plain way's flux is sampled at 8 phases
inside each bin here, so the bin comparison is at the level of that sampling (a few parts in 100,000); the comparison with
64 samples per bin is in likelihood_from_plain_star.py.
usage: python -I -B -X utf8 kernel_against_plain_star.py        (about six minutes; the apjlab5968 folder from settings.cfg)
The last line is RESULT ... (not part of check.py --slow; run on its own)."""
import math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk, plain_star as pst, quick_example as qe, settings
from consts import PC_CM
ref_dir = settings.ref_dir()
result = {}
for name, (ths, rho, inc, nu), ob, ga in (("SD1e", (60.0, 1.0, 30.0, 400.0), False, False), ("OS1e", (60.0, 1.0, 30.0, 600.0), True, True),
                                           ("OS1f", (20.0, 1.0, 80.0, 600.0), True, True)):
    d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
    y = pst.light_curve(1.4, 12.0, nu, 200.0 * PC_CM, math.radians(inc), math.radians(ths), rho, 0.35, [1.0], d[:, 0],
                        n_zeta=64, n_eta=256, n_phi=2048, oblate=ob, gamma_area=ga)[:, 0]
    e = float(np.max(np.abs(y - d[:, 1])) / np.median(d[:, 1]))
    result["plain_vs_" + name] = e
    print("the plain way against the agreed reference curve %s: largest difference over the median %.1e" % (name, e), flush=True)
NB, SUB = 32, 8
edges = np.linspace(0.0, 1.0, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
E = np.array([0.2, 0.5, 1.0, 1.6])
P2 = dict(qe.P, radius=10.0, log10T=6.8)
for key, label, P in (("quick", "the quick example's star", qe.P), ("second", "the second check star", P2)):
    M, D, kT, inc = qe.star_args(P)
    kern = {}
    for nr, n in ((10, 16), (24, 24), (48, 32)):
        kern[(nr, n)] = fk.spot_counts(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, edges, phase0=0.0, n_rings=nr, n=n)
    for nz, ne, nphi in ((48, 192, 2048), (96, 384, 4096)):
        t0 = time.perf_counter()
        y = pst.light_curve(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, mid, n_zeta=nz, n_eta=ne, n_phi=nphi)
        yb = y.reshape(NB, SUB, E.size).mean(axis=1) / NB
        # the sampled mean inside a bin is not the bin's integral exactly, so the differences here include the sampling
        for k in ((10, 16), (24, 24), (48, 32)):
            d = (kern[k] - yb) / yb.mean(axis=0)[None, :]
            print("%s: kernel (%d rings per panel, %d points) minus the plain way (%d x %d patches): largest %.1e, r.m.s. %.1e of the mean bin ; per energy largest %s  (%.0f s)"
                  % (label, k[0], k[1], nz, ne, np.max(np.abs(d)), np.sqrt(np.mean(d**2)), np.array2string(np.max(np.abs(d), axis=0), precision=1), time.perf_counter() - t0), flush=True)
            if k == (48, 32):
                result["%s_%dx%d_bins" % (key, nz, ne)] = float(np.max(np.abs(d)))
print("RESULT " + " ".join("%s=%.3e" % kv for kv in sorted(result.items())))
