"""Stars drawn at random from the public quick example's prior (mass 1.0 to 1.6, radius 10 to 13, distance 0.5 to 2,
cos i 0 to 1, phase shift -0.25 to 0.75, spot colatitude and radius 0.001 to pi/2 - 0.001, log10 T 6 to 7), worked out by
the Python program and by the compiled copy: does either stop, do the bins agree, do the likelihoods agree."""
import math, os, subprocess, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk, quick_example as qe, public_likelihood as pl
from consts import BB_CONST
counts, bg, support = qe.load_data(sys.argv[1])
ndraw = int(sys.argv[2]) if len(sys.argv) > 2 else 40
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph = np.linspace(0.0, 1.0, qe.NB + 1)
cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")
wsl_c = (lambda p: "/mnt/" + p[0].lower() + p[2:].replace(chr(92), "/") if p[1:2] == ":" else p)(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")))  # the same folder as the Linux side of a Windows machine sees it
npp, n, nE = 10, 16, 24
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
rng = np.random.default_rng(77)
worst_bin, worst_l, n_py_stop, n_c_stop, n_finite = 0.0, 0.0, 0, 0, 0
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
    dump = "dump_prior_%d.bin" % k
    cmd = "cd %s && ./nicer_fast tables.bin %d %d %d 0 %s %s" % (wsl_c, npp, n, nE, " ".join(repr(float(x)) for x in p), dump)
    out = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", cmd], capture_output=True, text=True).stdout.replace(chr(0), "")
    line = [l for l in out.splitlines() if l.startswith("ln L")]
    if not line:
        n_c_stop += 1
        print("draw %d: the compiled copy said: %s ; Python ln L %s ; star %s" % (k, out.strip()[:60], v, np.round(p, 4).tolist()))
        continue
    vc = float(line[0].split("=")[1])
    raw = np.fromfile(os.path.join(cdir, dump), dtype=np.float64)
    os.remove(os.path.join(cdir, dump))
    Sc = raw[:qe.NB * nE].reshape(qe.NB, nE)
    eS = np.max(np.abs(Sc - S)) / max(S.mean(), 1e-300)
    worst_bin = max(worst_bin, eS)
    if np.isfinite(v) and np.isfinite(vc):
        n_finite += 1
        worst_l = max(worst_l, abs(vc - v) / max(1.0, abs(v)))
    elif np.isfinite(v) != np.isfinite(vc):
        print("draw %d: one program left by the early exit and the other did not: Python %s, compiled %s" % (k, v, vc))
    if eS > 1e-7:
        print("draw %d: bins differ by %.1e of the mean bin ; star %s" % (k, eS, np.round(p, 4).tolist()))
print("%d draws from the prior: the Python program stopped on %d, the compiled copy on %d ; largest bin difference %.1e of the mean bin ; %d draws with a finite likelihood in both, largest relative difference %.1e"
      % (ndraw, n_py_stop, n_c_stop, worst_bin, n_finite, worst_l))
