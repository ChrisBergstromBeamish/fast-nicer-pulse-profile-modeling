# -*- coding: utf-8 -*-
"""Our programs beside the public code's measured times and values (the .jsonl file written by time_public_code.py inside
X-PSI's environment: run_public_code_timing.sh writes results/public_code_times.jsonl and refine_public_code.sh writes
results/public_code_refine.jsonl, both in 2-python-solution). For every setting of the public code: its value at the
check's star; where its bins were saved beside the .jsonl file, the bins side by side and its count rates put through our
likelihood; ours converged at the same points near the star, and the differences; and, for the quick check only, the
public code's time, ours in Python and compiled, and the ratio of the times.

Two checks are handled: "quick" (the shipped quick example: background limited to 0.9 to 1.1 of a given one) and "second"
(xpsi/tests/test_likelihood.py: radius 10, log10 T 6.8, no upper limit on the background).

Changed on 23 Nov 2026, after the first run of the public code on this machine:
  - the public code's points that fell outside its own limits are left out, and the printout says how many. At the
    second check the star sits on the public code's lower limit of the radius, so about half of the points near it that
    time_public_code.py picks are refused at once (the public code answers with no evaluation); their times are not times
    of an evaluation;
  - times and ratios are given for the quick check only; for the second check, values and bins only;
  - ours in Python is timed as the public code was: one evaluation per point, a different star each time, over the same
    points (five passes), fastest and median; the compiled program is timed by its own mode -300 (a different star at
    every working-out, the mean over 300 of them, fastest of 7 rounds).

usage: python -I -B -X utf8 2-python-solution/tools/compare_with_public_code.py <public_code_times.jsonl> [a note on the power state]
       The public example's data folder comes from settings.cfg; the compiled rows need `make` run in 3-compiled-solution.
       Run it when nothing else is running on the machine, and say in the note whether the machine is on battery or on
       mains power (on this laptop every time was longer on battery).
"""
import json
import math
import os
import sys
import time
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_v] = "1"          # one thread, set before NumPy loads, so that the time is for one processor core
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fast_kernel as fk
import public_likelihood as pl
import quick_example as qe
import settings
from consts import BB_CONST

if len(sys.argv) < 2:
    print("usage: python -I -B -X utf8 2-python-solution/tools/compare_with_public_code.py <public_code_times.jsonl> [a note on the power state]")
    sys.exit(2)
counts, bg, support = qe.load_data(settings.data_dir())
nosup = np.stack([np.zeros(counts.shape[0]), -np.ones(counts.shape[0])], axis=1)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph_edges = np.linspace(0.0, 1.0, qe.NB + 1)


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


def compiled(p, npp, n, nE, reps=0):
    """Our compiled copy (3-compiled-solution/nicer_fast, run through settings.run_compiled), with the quick check's limits
    on the background: its ln L at the star p, and with reps its time per working-out (s)."""
    out = settings.run_compiled("nicer_fast", ["tables.bin", npp, n, nE, reps] + [repr(float(x)) for x in p])
    val, tm = None, None
    for line in out.splitlines():
        if line.startswith("ln L ="):
            val = float(line.split("=")[1])
        if line.startswith("time per working-out"):
            tm = float(line.split(":")[1].split("ms")[0]) * 1e-3
    return val, tm


def times_over_points(f, pts, passes=5):
    """One evaluation per point, a different star each time, as the public code was timed; all the single times."""
    f(pts[0])
    out = []
    for _ in range(passes):
        for q in pts:
            t0 = time.perf_counter()
            f(q)
            out.append(time.perf_counter() - t0)
    return sorted(out)


SETTINGS = (("10 rings per panel, 16 points, 24 energies", (10, 16, 24)), ("8 rings per panel, 12 points, 16 energies", (8, 12, 16)))
P_QUICK = [1.4, 12.0, 1.0, math.cos(60 * math.pi / 180), 0.0, 70 * math.pi / 180, 0.75, 6.7]
recs = [json.loads(line) for line in open(sys.argv[1], encoding="utf-8") if line.strip()]
if len(sys.argv) > 2:
    print("power state when this was run: %s" % sys.argv[2])
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
    # the bins, where the public code's were saved
    npz = "%s.%s_cells%d_leaves%d_rays%d_energies%d.npz" % (sys.argv[1], mode, s["sqrt_num_cells"], s["num_leaves"], s["num_rays"], s["num_energies"])
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
    ok = [i for i, val in enumerate(rec["values"]) if val > -1.0e10]
    if len(ok) < len(rec["values"]):
        print("   %d of its %d points near the star were outside its own limits (it answered at once, with no evaluation); left out"
              % (len(rec["values"]) - len(ok), len(rec["values"])))
    pts = [rec["points"][i] for i in ok]
    if not pts:
        continue
    oc = np.array([conv(q) for q in pts])
    xv = np.array([rec["values"][i] for i in ok])
    print("   at its %d points near the star: public minus ours converged: mean %+.4f, largest %.4f"
          % (len(pts), float(np.mean(xv - oc)), float(np.max(np.abs(xv - oc)))))
    if mode != "quick":
        continue                  # times and ratios for the quick check only (the docstring says why)
    xt = sorted(rec["times_s"][i] for i in ok)
    x_fast, x_med = xt[0], xt[len(xt) // 2]
    print("   its time per likelihood evaluation, a different star each time: fastest %.4f s, median %.4f s over %d evaluations"
          % (x_fast, x_med, len(xt)))
    for label, (npp, n, nE) in SETTINGS:
        f = make(npp, n, nE, sup)
        ow = np.array([f(q) for q in pts])
        ts = times_over_points(f, pts)
        t_fast, t_med = ts[0], ts[len(ts) // 2]
        print("   ours in Python, %s: largest difference from ours converged %.4f ; fastest %.2f ms, median %.2f ms over %d evaluations ; public code's time / ours: fastest %.0f, median %.0f"
              % (label, float(np.max(np.abs(ow - oc))), t_fast * 1e3, t_med * 1e3, len(ts), x_fast / t_fast, x_med / t_med))
    for label, (npp, n, nE) in SETTINGS:
        vals = [compiled(q, npp, n, nE)[0] for q in pts]
        if any(v is None for v in vals):
            print("   ours compiled, %s: no answer from the compiled program (is it built?)" % label)
            continue
        tm = compiled(pts[0], npp, n, nE, reps=-300)[1]
        print("   ours compiled, %s: largest difference from ours converged %.4f ; %.3f ms per evaluation (mean over 300 different stars, fastest of 7 rounds) ; public code's time / ours: its fastest %.0f, its median %.0f"
              % (label, float(np.max(np.abs(np.array(vals) - oc))), tm * 1e3, x_fast / tm, x_med / tm))
