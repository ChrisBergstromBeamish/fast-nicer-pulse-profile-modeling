# -*- coding: utf-8 -*-
"""Our programs beside the public code's measured times and values (results/public_code_times.jsonl, written by
time_public_code.py inside X-PSI's environment). For every setting of the public code: its value at the check's star and
its time; then ours at the same points: values, differences, times, and the ratio of the times; and, where the public
code's bins were saved, the bins side by side.

Two checks are handled: "quick" (the shipped quick example: background limited to 0.9 to 1.1 of a given one) and "second"
(xpsi/tests/test_likelihood.py: radius 10, log10 T 6.8, no upper limit on the background).

usage: python compare_with_public_code.py <the example's Data folder> <public_code_times.jsonl>
"""
import json
import math
import os
import subprocess
import sys
import time
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_v] = "1"          # one thread, set before NumPy loads, so that the time is for one processor core
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
import public_likelihood as pl
import quick_example as qe
from consts import BB_CONST

counts, bg, support = qe.load_data(sys.argv[1])
nosup = np.stack([np.zeros(counts.shape[0]), -np.ones(counts.shape[0])], axis=1)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph_edges = np.linspace(0.0, 1.0, qe.NB + 1)
WSL_C = (lambda p: "/mnt/" + p[0].lower() + p[2:].replace(chr(92), "/") if p[1:2] == ":" else p)(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")))  # the same folder as the Linux side of a Windows machine sees it


def make(nr, n, nE, sup, want_star=False):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W

    def like(p):
        mass, radius, distance, cos_i, shift, colat, rho, logT = p
        M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
        S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                           qe.X_KB_OVER_KEV * 10.0**logT, E, ph_edges, phase0=0.5 - shift, n_rings=nr, n=n)
        star = Wn @ S.T
        return star if want_star else pl.marginal_loglike_fast(counts, star, qe.T_EXP, sup, fact)
    return like


def compiled(p, npp, n, nE, reps=0, no_upper=False):
    """Our compiled copy (c/nicer_fast, on the Linux side): its ln L at the star p, and with reps its time per working-out (s)."""
    cmd = "cd %s && %s./nicer_fast tables.bin %d %d %d %d %s" % (WSL_C, "NICER_NO_UPPER_LIMIT=1 " if no_upper else "", npp, n, nE, reps,
                                                               " ".join(repr(float(x)) for x in p))
    out = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", cmd], capture_output=True, text=True).stdout.replace(chr(0), "")
    val, tm = None, None
    for line in out.splitlines():
        if line.startswith("ln L ="):
            val = float(line.split("=")[1])
        if line.startswith("time per working-out"):
            tm = float(line.split(":")[1].split("ms")[0]) * 1e-3
    return val, tm


def best_time(f, p, reps=5, inner=20):
    best = 1e9
    for _ in range(reps):
        t0 = time.perf_counter()
        for _ in range(inner):
            f(p)
        best = min(best, (time.perf_counter() - t0) / inner)
    return best


SETTINGS = (("10 rings per panel, 16 points, 24 energies", (10, 16, 24)), ("8 rings per panel, 12 points, 16 energies", (8, 12, 16)))
P_QUICK = [1.4, 12.0, 1.0, math.cos(60 * math.pi / 180), 0.0, 70 * math.pi / 180, 0.75, 6.7]
recs = [json.loads(line) for line in open(sys.argv[2], encoding="utf-8") if line.strip()]
conv_at = {}
for rec in recs:
    mode = rec.get("mode", "quick")
    sup = nosup if mode == "second" else support
    p_star = rec.get("star", P_QUICK)
    conv = make(48, 32, 48, sup)
    if mode not in conv_at:
        conv_at[mode] = conv(p_star)
        print("\n=== the %s check: ours, converged, at its star: %.6f ; the public code's shipped value %.6f ; ours minus shipped %+.6f"
              % (mode, conv_at[mode], rec["shipped"], conv_at[mode] - rec["shipped"]))
    s = rec["settings"]
    print("\npublic code at sqrt_num_cells %d, leaves %d, rays %d, energies %d (version %s), %s check:"
          % (s["sqrt_num_cells"], s["num_leaves"], s["num_rays"], s["num_energies"], rec.get("xpsi_version", "?"), mode))
    print("   at the star: %.6f ; minus its shipped value %+.6f ; minus ours converged %+.6f"
          % (rec["value_at_p"], rec["value_at_p"] - rec["shipped"], rec["value_at_p"] - conv_at[mode]))
    print("   time per likelihood evaluation: fastest %.4f s, median %.4f s over %d calls"
          % (rec["fastest_s"], rec["median_s"], len(rec["times_s"])))
    # the bins, where the public code's were saved
    npz = "%s.%s_cells%d_leaves%d_rays%d_energies%d.npz" % (sys.argv[2], mode, s["sqrt_num_cells"], s["num_leaves"], s["num_rays"], s["num_energies"])
    if os.path.isfile(npz):
        z = np.load(npz)
        if "rate" in z and "signal_phases" in z:
            ours_star = make(48, 32, 48, sup, want_star=True)(p_star)                      # mean count rate per channel and bin
            theirs = pl.phase_integrate(np.asarray(z["rate"]), np.asarray(z["signal_phases"]), ph_edges, p_star[4])
            rel = theirs / ours_star - 1.0
            print("   bins (its rates through its own phase step, against ours converged): total %+.2e ; largest %.2e ; r.m.s. %.2e ; per phase bin largest %.2e ; per channel largest %.2e"
                  % (theirs.sum() / ours_star.sum() - 1.0, np.max(np.abs(rel)), np.sqrt(np.mean(rel**2)),
                     np.max(np.abs(theirs.sum(axis=0) / ours_star.sum(axis=0) - 1.0)), np.max(np.abs(theirs.sum(axis=1) / ours_star.sum(axis=1) - 1.0))))
            # the exchange: its rates through our likelihood
            print("   its rates through our likelihood: %.6f (its own likelihood gave %.6f)"
                  % (pl.marginal_loglike_fast(counts, theirs, qe.T_EXP, sup, fact), rec["value_at_p"]))
    pts = rec["points"]
    if not pts:
        continue
    oc = np.array([conv(q) for q in pts])
    xv = np.array(rec["values"])
    print("   at its %d points near the star: public minus ours converged: mean %+.4f, largest %.4f"
          % (len(pts), float(np.mean(xv - oc)), float(np.max(np.abs(xv - oc)))))
    for label, (npp, n, nE) in SETTINGS:
        f = make(npp, n, nE, sup)
        ow = np.array([f(q) for q in pts])
        t = best_time(f, pts[0])
        print("   ours in Python, %s: largest difference from ours converged %.4f ; %.2f ms per evaluation ; public code's fastest time / ours = %.0f"
              % (label, float(np.max(np.abs(ow - oc))), t * 1e3, rec["fastest_s"] / t))
    for label, (npp, n, nE) in SETTINGS:
        vals = [compiled(q, npp, n, nE, no_upper=(mode == "second"))[0] for q in pts]
        if any(v is None for v in vals):
            print("   ours compiled, %s: no answer from the compiled program (is it built?)" % label)
            continue
        tm = compiled(pts[0], npp, n, nE, reps=-300, no_upper=(mode == "second"))[1]
        print("   ours compiled, %s: largest difference from ours converged %.4f ; %.3f ms per evaluation ; public code's fastest time / ours = %.0f"
              % (label, float(np.max(np.abs(np.array(vals) - oc))), tm * 1e3, rec["fastest_s"] / tm))
