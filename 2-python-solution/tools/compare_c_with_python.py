"""The compiled program (c/nicer_fast) against our Python programs: the star's bin counts S, the model rates per channel
and bin, and the likelihood, for a list of stars. The compiled program must have been run with a dump file for each star."""
import math, os, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk, quick_example as qe, public_likelihood as pl
from consts import BB_CONST
data_dir = sys.argv[1]
counts, bg, support = qe.load_data(data_dir)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph = np.linspace(0.0, 1.0, qe.NB + 1)
here = os.path.dirname(os.path.abspath(__file__))
cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")
wsl_c = (lambda p: "/mnt/" + p[0].lower() + p[2:].replace(chr(92), "/") if p[1:2] == ":" else p)(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")))  # the same folder as the Linux side of a Windows machine sees it
stars = [(1.4, 12.0, 1.0, math.cos(math.radians(60.0)), 0.0, math.radians(70.0), 0.75, 6.7),
         (1.4011, 12.007, 0.9995, 0.5003, 0.0007, 1.2224, 0.7494, 6.7003),
         (1.2, 12.8, 1.1, 0.8, 0.1, 0.5, 1.2, 6.69),
         (1.55, 10.4, 0.9, 0.1, -0.2, 1.5, 0.3, 6.72),
         (1.1, 11.0, 1.0, 0.999999, 0.3, 0.9, 0.6, 6.7),
         (1.3, 12.5, 1.0, 0.3, 0.6, 0.2, 0.9, 6.68)]
worstS, worstL = 0.0, 0.0
for (npp, n, nE) in ((10, 16, 24), (8, 12, 16), (24, 24, 32)):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W
    for k, p in enumerate(stars):
        mass, radius, distance, cos_i, shift, colat, rho, logT = p
        M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
        S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                           qe.X_KB_OVER_KEV * 10.0**logT, E, ph, phase0=0.5 - shift, n_rings=npp, n=n) * (qe.X_BB / BB_CONST)
        star = (qe.AREA * qe.NB) * (W @ S.T)
        v = pl.marginal_loglike_fast(counts, star, qe.T_EXP, support, fact)
        dump = "dump_%d_%d_%d_%d.bin" % (npp, n, nE, k)
        cmd = "cd %s && ./nicer_fast tables.bin %d %d %d 0 %s %s" % (wsl_c, npp, n, nE, " ".join(repr(float(x)) for x in p), dump)
        out = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", cmd], capture_output=True, text=True).stdout.replace("\x00", "")
        line = [l for l in out.splitlines() if l.startswith("ln L")]
        if not line:
            print("star %d at (%d, %d, %d): the compiled program said: %s ; Python ln L %s" % (k, npp, n, nE, out.strip()[:80], v))
            continue
        vc = float(line[0].split("=")[1])
        raw = np.fromfile(os.path.join(cdir, dump), dtype=np.float64)
        Sc = raw[:qe.NB * nE].reshape(qe.NB, nE)
        starc = raw[qe.NB * nE:].reshape(291, qe.NB)
        eS = np.max(np.abs(Sc - S)) / S.mean()
        eR = np.max(np.abs(starc - star)) / star.mean()
        dL = vc - v if np.isfinite(v) else (0.0 if not np.isfinite(vc) else float("nan"))
        worstS, worstL = max(worstS, eS), max(worstL, abs(dL))
        print("star %d at (%d rings per panel, %d points, %d energies): bin counts differ by %.1e of the mean bin, model rates by %.1e ; ln L compiled %.9f, Python %.9f, difference %+.1e"
              % (k, npp, n, nE, eS, eR, vc, v, dL), flush=True)
        os.remove(os.path.join(cdir, dump))
print("largest differences over all: bin counts %.1e of the mean bin ; ln L %.1e" % (worstS, worstL))
