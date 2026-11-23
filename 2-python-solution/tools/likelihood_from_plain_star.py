"""The likelihood of the public quick example and of its second check, worked with the plain way for flattened stars
(plain_star.py: patches, direction search, lensing factor written out; no count at emission, no rings, no function K), and
set beside the fast kernel's and the public code's shipped numbers. The plain way's flux is sampled at 64 phases inside each
bin; its accuracy is a few parts in 100,000 per bin, so at the second star (where 1 part in 10,000 moves ln L by 9.4) it can
only say which of the two numbers it is nearer to."""
import math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk, plain_star as pst, quick_example as qe, public_likelihood as pl
from consts import BB_CONST
counts, bg, support = qe.load_data(sys.argv[1])
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
nosup = np.stack([np.zeros(291), -np.ones(291)], axis=1)
NB, SUB, nE = 32, 64, 24
ph_edges = np.linspace(0.0, 1.0, NB + 1)
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W
P2 = dict(qe.P, radius=10.0, log10T=6.8)
for label, P, sup, shipped, tol in (("quick example's star", qe.P, support, -32552.616905, 0.3255), ("second check star", P2, nosup, -47881.27817666349, 0.479)):
    M, D, kT, inc = qe.star_args(P)
    Sk = fk.spot_counts(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, ph_edges, phase0=0.5, n_rings=48, n=32)
    vk = pl.marginal_loglike_fast(counts, Wn @ Sk.T, qe.T_EXP, sup, fact)
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
