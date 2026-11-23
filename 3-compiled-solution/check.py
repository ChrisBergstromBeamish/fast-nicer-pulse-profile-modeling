# -*- coding: utf-8 -*-
"""One command that compares the compiled solution, in both of its built forms, with the Python solution.

usage, from the repository's top, on the side that has Python with NumPy:
    python -I -B -X utf8 3-compiled-solution/check.py [nicer_fast | nicer_fast_strict]
Without a name both forms are checked (make in 3-compiled-solution builds them; on Windows they run on the Linux side through
WSL, see settings.cfg). For each form, four comparisons with the Python solution, eight lines:
  1. the public quick example at the working setting (10 rings per panel, 16 points, 24 energies): the likelihood
  2. six stars at three settings: the phase-bin counts, bin by bin, the likelihood, and that both programs take the public
     code's early exit on the same stars                                                    (tools/compare_c_with_python.py)
  3. sixty stars drawn at random from the public example's prior: the same, and that neither program stops
                                                                                            (tools/prior_draws_c_and_python.py)
  4. the likelihood step alone on 247 hard made-up channels, channel by channel             (tools/hard_channels_through_c.py)
Each line says: ok or FAIL, what is compared, the difference found, the limit. The exit code is 0 if every check holds, 1 if any
fails, 2 if a program is not built or the settings are missing. Every number checked here has a row in
1-setup-and-documentation/CLAIMS.md. The Python solution is the reference; this command does not time anything.
"""
import math
import os
import re
import subprocess
import sys
import time
HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "..", "2-python-solution", "tools")
sys.path.insert(0, os.path.join(HERE, "..", "2-python-solution", "code"))
sys.path.insert(0, TOOLS)
import settings

PROGRAMS = ("nicer_fast", "nicer_fast_strict")
wanted = sys.argv[1:] or list(PROGRAMS)
if any(p not in PROGRAMS for p in wanted):
    print("usage: python -I -B -X utf8 3-compiled-solution/check.py [nicer_fast | nicer_fast_strict]")
    sys.exit(2)
try:
    data_dir = settings.data_dir()
except settings.SettingsError as exc:
    print("cannot run the checks: %s" % exc)
    sys.exit(2)
import numpy as np                       # noqa: E402
import fast_kernel as fk                 # noqa: E402
import public_likelihood as pl           # noqa: E402
import quick_example as qe               # noqa: E402
from consts import BB_CONST              # noqa: E402

results = []
T0 = time.perf_counter()
print("Checks of the compiled solution against the Python solution. The public example's data: %s" % data_dir, flush=True)


def check(name, value, limit):
    ok = bool(value <= limit)
    results.append(ok)
    print("%s  %-110s %.2e (limit %.3g)" % ("ok  " if ok else "FAIL", name, value, limit), flush=True)


def script_result(script, args):
    """Run one tools/ script as a command, show its output indented, and read its 'RESULT key=value ...' line."""
    t1 = time.perf_counter()
    done = subprocess.run([sys.executable, "-I", "-B", "-X", "utf8", os.path.join(TOOLS, script)] + args, capture_output=True, text=True)
    out = (done.stdout or "") + (done.stderr or "")
    for line in out.splitlines():
        print("      " + line)
    print("      [%s: %.0f s]" % (script, time.perf_counter() - t1), flush=True)
    found = {}
    for line in out.splitlines():
        if line.startswith("RESULT "):
            for key, value in re.findall(r"(\w+)=(\S+)", line):
                found[key] = float(value)
    return found


# the Python solution's likelihood of the quick example at the working setting (as 2-python-solution/check.py computes it)
counts, bg, support = qe.load_data(data_dir)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], 24, edges_E)
M, D, kT, inc = qe.star_args(qe.P)
S = fk.spot_counts(M, qe.P["radius"], qe.P["freq"], D, inc, qe.P["colat"], qe.P["ang_radius"], kT, E, np.linspace(0.0, 1.0, qe.NB + 1),
                   phase0=0.5 - qe.P["phase_shift"], n_rings=10, n=16)
v_py = pl.marginal_loglike_fast(counts, (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * (W @ S.T), qe.T_EXP, support, fact)
missing = []
for program in wanted:
    print("\n--- %s (%s)" % (program, "the compiler's fast arithmetic" if program == "nicer_fast" else "strict arithmetic"), flush=True)
    P = qe.P
    out = settings.run_compiled(program, ["tables.bin", 10, 16, 24, 0] + [repr(float(P[k])) for k in ("mass", "radius", "distance", "cos_i", "phase_shift", "colat", "ang_radius", "log10T")])
    line = [l for l in out.splitlines() if l.startswith("ln L =")]
    if not line:
        print("      %s gave no answer: is it built (make in 3-compiled-solution), is tables.bin there, is settings.cfg right? It said: %s"
              % (program, out.strip()[:200] or "(nothing)"))
        missing.append(program)
        continue
    check("%s: the quick example at the working setting, ln L against the Python solution" % program, abs(float(line[0].split("=")[1]) - v_py), 1e-6)
    r = script_result("compare_c_with_python.py", [program])
    check("%s: six stars at three settings, phase-bin counts, of the mean bin (one star seen 0.08 degrees from its axis; the rest 1e-10)" % program,
          r.get("largest_bins", math.inf), 1e-7)
    check("%s: six stars at three settings, ln L, relative (|difference| / |ln L|)" % program, r.get("largest_lnL_rel", math.inf), 1e-8)
    check("%s: six stars at three settings, cases where only one program took the public code's early exit (ln L = -inf)" % program, r.get("mismatch", math.inf), 0)
    r = script_result("prior_draws_c_and_python.py", [program, "60"])
    check("%s: sixty prior draws, workings-out on which either program stopped, or only one took the early exit" % program,
          r.get("py_stopped", math.inf) + r.get("c_stopped", math.inf) + r.get("mismatch", math.inf), 0)
    check("%s: sixty prior draws, phase-bin counts, of the mean bin" % program, r.get("largest_bins", math.inf), 1e-7)
    check("%s: sixty prior draws with a finite likelihood in both, ln L, relative" % program, r.get("largest_lnL_rel", math.inf), 1e-9)
    r = script_result("hard_channels_through_c.py", [program])
    check("%s: the likelihood step on 247 hard made-up channels against the step-by-step Python form, per channel" % program, r.get("largest", math.inf), 1e-8)
print("\n%d checks: %d hold, %d fail%s ; %.0f s" % (len(results), sum(results), len(results) - sum(results),
      " ; not built or not run: %s" % ", ".join(missing) if missing else "", time.perf_counter() - T0))
sys.exit(2 if missing else (0 if all(results) else 1))
