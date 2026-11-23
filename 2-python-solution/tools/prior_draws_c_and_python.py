"""Stars drawn at random from the public quick example's prior (mass 1.0 to 1.6, radius 10 to 13, distance 0.5 to 2,
cos i 0 to 1, phase shift -0.25 to 0.75, spot colatitude and radius 0.001 to pi/2 - 0.001, log10 T 6 to 7), worked out by
the Python program and by the compiled copy: does either stop, do the bins agree, do the likelihoods agree.
usage: python -I -B -X utf8 prior_draws_c_and_python.py [nicer_fast | nicer_fast_strict] [number of draws]   (nicer_fast, 40)
The last line, RESULT ..., is read by 3-compiled-solution/check.py."""
import math, os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk, quick_example as qe, public_likelihood as pl, settings
from consts import BB_CONST
program = sys.argv[1] if len(sys.argv) > 1 else "nicer_fast"
ndraw = int(sys.argv[2]) if len(sys.argv) > 2 else 40
counts, bg, support = qe.load_data(settings.data_dir())
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph = np.linspace(0.0, 1.0, qe.NB + 1)
cdir = settings.compiled_dir()
npp, n, nE = 10, 16, 24
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
rng = np.random.default_rng(77)
worst_bin, worst_l, n_py_stop, n_c_stop, n_finite, n_mismatch = 0.0, 0.0, 0, 0, 0, 0
for k in range(ndraw):
    p = [rng.uniform(1.0, 1.6), rng.uniform(10.0, 13.0), rng.uniform(0.5, 2.0), rng.uniform(0.0, 1.0), rng.uniform(-0.25, 0.75),
         rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(6.0, 7.0)]
    mass, radius, distance, cos_i, shift, colat, rho, logT = p
    M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
    S = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                               qe.X_KB_OVER_KEV * 10.0**logT, E, ph, phase0=0.5 - shift, n_rings=npp, n=n) * (qe.X_BB / BB_CONST)
        v = pl.marginal_loglike_fast(counts, (qe.AREA * qe.NB) * (W @ S.T), qe.T_EXP, support, fact)
    except Exception as ex:
        n_py_stop += 1
        print("draw %d: the Python program stopped: %s: %s ; star %s" % (k, type(ex).__name__, str(ex)[:80], np.round(p, 4).tolist()))
        continue
    dump = "dump_prior_%d_%d.bin" % (k, os.getpid())        # the process number keeps two runs apart
    out = settings.run_compiled(program, ["tables.bin", npp, n, nE, 0] + [repr(float(x)) for x in p] + [dump])
    line = [l for l in out.splitlines() if l.startswith("ln L")]
    if not line:
        n_c_stop += 1
        print("draw %d: the compiled copy said: %s ; Python ln L %s ; star %s" % (k, out.strip()[:60], v, np.round(p, 4).tolist()))
        continue
    vc = line[0].split("=")[1].strip()
    vc = -math.inf if vc == "-inf" else float(vc)
    raw = np.fromfile(os.path.join(cdir, dump), dtype=np.float64)
    os.remove(os.path.join(cdir, dump))
    Sc = raw[:qe.NB * nE].reshape(qe.NB, nE)
    eS = np.max(np.abs(Sc - S)) / max(S.mean(), 1e-300)
    worst_bin = max(worst_bin, eS)
    if np.isfinite(v) and np.isfinite(vc):
        n_finite += 1
        worst_l = max(worst_l, abs(vc - v) / max(1.0, abs(v)))
    elif np.isfinite(v) != np.isfinite(vc):
        n_mismatch += 1
        print("draw %d: one program left by the early exit and the other did not: Python %s, compiled %s" % (k, v, vc))
    if eS > 1e-7:
        print("draw %d: bins differ by %.1e of the mean bin ; star %s" % (k, eS, np.round(p, 4).tolist()))
print("%d draws from the prior: the Python program stopped on %d, the compiled copy on %d ; largest bin difference %.1e of the mean bin ; %d draws with a finite likelihood in both, largest relative difference %.1e ; early exit in one program only: %d"
      % (ndraw, n_py_stop, n_c_stop, worst_bin, n_finite, worst_l, n_mismatch))
print("RESULT draws=%d py_stopped=%d c_stopped=%d mismatch=%d largest_bins=%.3e finite=%d largest_lnL_rel=%.3e"
      % (ndraw, n_py_stop, n_c_stop, n_mismatch, worst_bin, n_finite, worst_l))
