"""The compiled program (3-compiled-solution/nicer_fast, or nicer_fast_strict) against our Python programs: the star's bin
counts S, the model rates per channel and bin, and the likelihood, for six stars at three settings. The compiled program writes
a dump file for each star, which is read back and removed. It runs on the Linux side (settings.py).
usage: python -I -B -X utf8 compare_c_with_python.py [nicer_fast | nicer_fast_strict]     (nicer_fast if no name is given)
The last line, RESULT ..., is read by 3-compiled-solution/check.py."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk, quick_example as qe, public_likelihood as pl, settings
from consts import BB_CONST
program = sys.argv[1] if len(sys.argv) > 1 else "nicer_fast"
data_dir = settings.data_dir()
counts, bg, support = qe.load_data(data_dir)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph = np.linspace(0.0, 1.0, qe.NB + 1)
cdir = settings.compiled_dir()
stars = [(1.4, 12.0, 1.0, math.cos(math.radians(60.0)), 0.0, math.radians(70.0), 0.75, 6.7),
         (1.4011, 12.007, 0.9995, 0.5003, 0.0007, 1.2224, 0.7494, 6.7003),
         (1.2, 12.8, 1.1, 0.8, 0.1, 0.5, 1.2, 6.69),
         (1.55, 10.4, 0.9, 0.1, -0.2, 1.5, 0.3, 6.72),
         (1.1, 11.0, 1.0, 0.999999, 0.3, 0.9, 0.6, 6.7),
         (1.3, 12.5, 1.0, 0.3, 0.6, 0.2, 0.9, 6.68)]
worstS, worstR, worstL, worstLrel, n_missing, n_mismatch = 0.0, 0.0, 0.0, 0.0, 0, 0
for (npp, n, nE) in ((10, 16, 24), (8, 12, 16), (24, 24, 32)):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    for k, p in enumerate(stars):
        mass, radius, distance, cos_i, shift, colat, rho, logT = p
        M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
        S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                           qe.X_KB_OVER_KEV * 10.0**logT, E, ph, phase0=0.5 - shift, n_rings=npp, n=n) * (qe.X_BB / BB_CONST)
        star = (qe.AREA * qe.NB) * (W @ S.T)
        v = pl.marginal_loglike_fast(counts, star, qe.T_EXP, support, fact)
        dump = "dump_%d_%d_%d_%d_%d.bin" % (npp, n, nE, k, os.getpid())        # the process number keeps two runs apart
        out = settings.run_compiled(program, ["tables.bin", npp, n, nE, 0] + [repr(float(x)) for x in p] + [dump])
        line = [l for l in out.splitlines() if l.startswith("ln L")]
        if not line:
            n_missing += 1
            print("star %d at (%d, %d, %d): the compiled program said: %s ; Python ln L %s" % (k, npp, n, nE, out.strip()[:80], v))
            continue
        vc = line[0].split("=")[1].strip()
        vc = -math.inf if vc == "-inf" else float(vc)
        raw = np.fromfile(os.path.join(cdir, dump), dtype=np.float64)
        Sc = raw[:qe.NB * nE].reshape(qe.NB, nE)
        starc = raw[qe.NB * nE:].reshape(291, qe.NB)
        eS = np.max(np.abs(Sc - S)) / S.mean()
        eR = np.max(np.abs(starc - star)) / star.mean()
        if np.isfinite(v) and np.isfinite(vc):
            dL = vc - v
            rel = abs(dL) / max(1.0, abs(v))
        elif np.isfinite(v) or np.isfinite(vc):                # the public code's early exit in one program only
            n_mismatch += 1
            dL = rel = math.inf
        else:
            dL = rel = 0.0                                     # both left by the early exit
        worstS, worstR, worstL, worstLrel = max(worstS, eS), max(worstR, eR), max(worstL, abs(dL)), max(worstLrel, rel)
        print("star %d at (%d rings per panel, %d points, %d energies): bin counts differ by %.1e of the mean bin, model rates by %.1e ; ln L compiled %.9f, Python %.9f, difference %+.1e"
              % (k, npp, n, nE, eS, eR, vc, v, dL), flush=True)
        os.remove(os.path.join(cdir, dump))
print("largest differences over all: bin counts %.1e of the mean bin ; model rates %.1e ; ln L %.1e, relative %.1e ; stars without an answer from the compiled program: %d ; early exit in one program only: %d"
      % (worstS, worstR, worstL, worstLrel, n_missing, n_mismatch))
print("RESULT largest_bins=%.3e largest_rates=%.3e largest_lnL=%.3e largest_lnL_rel=%.3e missing=%d mismatch=%d"
      % (worstS, worstR, worstL if n_missing == 0 else math.inf, worstLrel if n_missing == 0 else math.inf, n_missing, n_mismatch))
