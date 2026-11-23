"""The fast kernel itself against the sixteen agreed reference curves of the test article (Bogdanov et al. 2019, its
supplementary files): for each reference phase the kernel's mean flux in a narrow bin around it (1/2048 of a turn) is set
beside the reference value. Largest difference over the curve's median, as the article measures it."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
from consts import PC_CM
ref_dir = sys.argv[1]
SD1 = {"SD1a": (90.0, 0.01, 90.0, 1.0), "SD1b": (90.0, 1.0, 90.0, 1.0), "SD1c": (90.0, 0.01, 90.0, 200.0),
       "SD1d": (90.0, 1.0, 90.0, 200.0), "SD1e": (60.0, 1.0, 30.0, 400.0), "SD1f": (20.0, 1.0, 80.0, 400.0)}
OS1 = {"OS1a": (90.0, 0.01, 90.0, 600.0, None), "OS1b": (90.0, 1.0, 90.0, 600.0, None), "OS1c": (90.0, 0.01, 90.0, 200.0, None),
       "OS1d": (90.0, 1.0, 90.0, 1.0, None), "OS1e": (60.0, 1.0, 30.0, 600.0, None), "OS1f": (20.0, 1.0, 80.0, 600.0, None),
       "OS1g": (60.0, 1.0, 30.0, 600.0, "c"), "OS1h": (60.0, 1.0, 30.0, 600.0, "s"),
       "OS1i": (20.0, 1.0, 80.0, 600.0, "c"), "OS1j": (20.0, 1.0, 80.0, 600.0, "s")}
BEAM = {None: None, "c": lambda m: m * m, "s": lambda m: 1.0 - m * m}
w = 1.0 / 2048.0
for fam, tests, ob, ga in (("round", SD1, False, False), ("flattened", OS1, True, True)):
    for name, tup in tests.items():
        ths, rho, inc, nu = tup[:4]
        bm = BEAM[tup[4]] if len(tup) > 4 else None
        d = np.loadtxt(os.path.join(ref_dir, name + "_test_IM.txt"))
        out = []
        for nr, n in ((10, 16), (24, 24)):
            y = np.empty(d.shape[0])
            # every reference phase gets its own narrow bin; done in blocks so that the edges increase
            for start in range(0, d.shape[0], 1):
                pass
            edges = np.empty(2 * d.shape[0])
            edges[0::2] = d[:, 0] - 0.5 * w
            edges[1::2] = d[:, 0] + 0.5 * w
            S = fk.spot_counts(1.4, 12.0, nu, 200.0 * PC_CM, math.radians(inc), math.radians(ths), rho, 0.35, [1.0], edges,
                               n_rings=nr, n=n, oblate=ob, gamma_area=ga, beaming=bm)[:, 0]
            y = S[0::2] / w
            out.append(np.max(np.abs(y - d[:, 1])) / np.median(d[:, 1]))
        print("%s (%s star, spot radius %.2f): kernel minus the agreed curve, largest over the median: %.1e at 10 rings per panel and 16 points ; %.1e at 24 and 24"
              % (name, fam, rho, out[0], out[1]), flush=True)
