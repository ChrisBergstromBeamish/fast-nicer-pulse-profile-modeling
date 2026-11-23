"""The likelihood of the public quick example and of its second check, worked with the plain way for flattened stars
(plain_star.py: patches, direction search, lensing factor written out; no count at emission, no rings, no function K), and
set beside the fast kernel's and the public code's shipped numbers. The plain way's flux is sampled at 64 phases inside each
bin; its accuracy is a few parts in 100,000 per bin at 48 x 192 patches and a few parts in a million at 96 x 384, so at the
second star (where 1 part in 10,000 moves ln L by 9.4) it can only say which of the two numbers it is nearer to.
usage: python -I -B -X utf8 likelihood_from_plain_star.py        (about eight minutes; the public example's data from settings.cfg)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
import math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk, plain_star as pst, quick_example as qe, public_likelihood as pl, settings
from consts import BB_CONST
counts, bg, support = qe.load_data(settings.data_dir())
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
nosup = np.stack([np.zeros(291), -np.ones(291)], axis=1)
NB, SUB, nE = 32, 64, 24
ph_edges = np.linspace(0.0, 1.0, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W
P2 = dict(qe.P, radius=10.0, log10T=6.8)
result = {}
for key, label, P, sup, shipped, tol in (("quick", "quick example's star", qe.P, support, -32552.616905, 0.3255),
                                         ("second", "second check star", P2, nosup, -47881.27817666349, 0.479)):
    M, D, kT, inc = qe.star_args(P)
    Sk = fk.spot_counts(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, ph_edges, phase0=0.5, n_rings=48, n=32)
    vk = pl.marginal_loglike_fast(counts, Wn @ Sk.T, qe.T_EXP, sup, fact)
    result[key + "_kernel_lnL"] = vk
    for nz, ne, nphi in ((48, 192, 2048), (96, 384, 4096)):
        t0 = time.perf_counter()
        y = pst.light_curve(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, mid - 0.5, n_zeta=nz, n_eta=ne, n_phi=nphi)
        Sp = y.reshape(NB, SUB, nE).mean(axis=1) / NB
        vp = pl.marginal_loglike_fast(counts, Wn @ Sp.T, qe.T_EXP, sup, fact)
        dbin = (Sk - Sp) / Sp.mean(axis=0)[None, :]
        print("%s, plain way with %d x %d patches: ln L = %.4f ; the fast kernel %.4f ; the public code's shipped %.4f (tolerance %.3f)"
              % (label, nz, ne, vp, vk, shipped, tol))
        print("     plain minus kernel %+.4f ; plain minus shipped %+.4f ; kernel minus shipped %+.4f ; bins: kernel minus plain largest %.1e, r.m.s. %.1e of the mean bin  (%.0f s)"
              % (vp - vk, vp - shipped, vk - shipped, np.max(np.abs(dbin)), np.sqrt(np.mean(dbin**2)), time.perf_counter() - t0), flush=True)
        tag = "%s_%dx%d" % (key, nz, ne)
        result[tag + "_bins"] = float(np.max(np.abs(dbin)))
        result[tag + "_lnL"] = vp
        result[tag + "_minus_kernel_abs"] = abs(vp - vk)
print("RESULT " + " ".join("%s=%.4e" % (k, v) if k.endswith(("_bins", "_abs")) else "%s=%.4f" % (k, v) for k, v in sorted(result.items())))
