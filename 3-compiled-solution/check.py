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
Then the batch mode (many stars at once), checked against the program itself, not against Python (batch_check.sh):
  5. for each form (and nicer_fast_no_openmp, the fast form built without OpenMP), the batch gives every star exactly the bits it
     gets alone (ln L, the code, every bin count and model rate): the sixty prior draws, the six stars at three settings, and
     refused stars among others at three settings, on 1, 2, 3, 4 and 8 threads, and again with the stars in a shuffled order
  6. the same batch run twice gives the same bits
  7. the shared library through nicer_fast_lib.py (Python's ctypes, on the Linux side) gives the batch mode's lines exactly
  8. the race detector (make race: 256 stars on 8 threads) reports no race
Each line says: ok or FAIL, what is compared, the difference found, the limit. The exit code is 0 if every check holds, 1 if any
fails, 2 if a program is not built or the settings are missing. Every number checked here has a row in
1-setup-and-documentation/CLAIMS.md. The Python solution is the reference; this command does not time anything.
"""
import math
import os
import re
import shutil
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


def linux(command, timeout=900):
    """Run a shell command in 3-compiled-solution on the Linux side (directly on Linux, through WSL on Windows, as
    settings.run_compiled does) and return its exit code and what it printed."""
    cdir = settings.compiled_dir()
    if os.name == "posix":
        done = subprocess.run(["bash", "-c", command], cwd=cdir, capture_output=True, text=True, timeout=timeout)
    else:
        distro = settings.read_settings().get("wsl", "")
        argv = ["wsl.exe"] + (["-d", distro] if distro else []) + ["--", "bash", "-lc",
                                                                   "cd %s && %s" % (settings.linux_side_path(cdir), command)]
        done = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    return done.returncode, ((done.stdout or "") + (done.stderr or "")).replace("\x00", "")


def read_lines(path):
    with open(path, encoding="ascii") as fh:
        return [l.split() for l in fh.read().splitlines() if l.strip()]


# 5 and 6: the batch mode against each star alone, bit for bit (stars files written here, worked out by batch_check.sh)
cdir = settings.compiled_dir()
rng = np.random.default_rng(77)                                     # the sixty stars of tools/prior_draws_c_and_python.py
sets = {"sixty": [[rng.uniform(1.0, 1.6), rng.uniform(10.0, 13.0), rng.uniform(0.5, 2.0), rng.uniform(0.0, 1.0), rng.uniform(-0.25, 0.75),
                   rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(6.0, 7.0)] for k in range(60)],
        "six": [(1.4, 12.0, 1.0, math.cos(math.radians(60.0)), 0.0, math.radians(70.0), 0.75, 6.7),     # those of tools/compare_c_with_python.py
                (1.4011, 12.007, 0.9995, 0.5003, 0.0007, 1.2224, 0.7494, 6.7003), (1.2, 12.8, 1.1, 0.8, 0.1, 0.5, 1.2, 6.69),
                (1.55, 10.4, 0.9, 0.1, -0.2, 1.5, 0.3, 6.72), (1.1, 11.0, 1.0, 0.999999, 0.3, 0.9, 0.6, 6.7), (1.3, 12.5, 1.0, 0.3, 0.6, 0.2, 0.9, 6.68)],
        # refused stars (code -1: too compact, too large, spot radius 0, spot on the pole, distance 0) among stars that are
        # worked out (the quick example, one of the six) and stars that take the public code's early exit (ln L = -inf)
        "refused": [(1.6, 6.0, 1.0, 0.5, 0.0, 1.0, 0.5, 6.5), (1.4, 12.0, 1.0, 0.5000000000000001, 0.0, 1.2217304763960306, 0.75, 6.7),
                    (1.4, 30.0, 1.0, 0.5, 0.0, 1.0, 0.5, 6.5), (1.6, 10.0, 0.5, 0.0, 0.0, 1.4, 1.5, 7.0), (1.4, 12.0, 1.0, 0.5, 0.0, 1.2, 0.0, 6.5),
                    (1.0, 4.5, 1.0, 0.0, 0.0, 1.57, 0.1, 6.5), (1.4, 12.0, 0.02, 0.5, 0.0, 1.2, 0.75, 7.0), (1.4, 12.0, 1.0, 0.5, 0.0, 0.0, 0.75, 6.5),
                    (1.2, 12.8, 1.1, 0.8, 0.1, 0.5, 1.2, 6.69), (1.4, 12.0, 0.0, 0.5, 0.0, 1.2, 0.75, 6.5)]}
order = {}
for name, stars in sets.items():
    order[name] = np.random.default_rng(5).permutation(len(stars))
    for suffix, rows in (("", stars), ("_shuffled", [stars[i] for i in order[name]])):
        with open(os.path.join(cdir, "stars_check_%s%s.txt" % (name, suffix)), "w", encoding="ascii", newline="\n") as fh:
            fh.writelines(" ".join(repr(float(x)) for x in p) + "\n" for p in rows)
jobs = [((10, 16, 24), "sixty")] + [(s, name) for s in ((10, 16, 24), (8, 12, 16), (24, 24, 32)) for name in ("six", "refused")]
forms = wanted + (["nicer_fast_no_openmp"] if "nicer_fast" in wanted else [])
for program in forms:
    print("\n--- %s: the batch mode against each star alone (batch_check.sh)" % program, flush=True)
    t1 = time.perf_counter()
    code, out = linux("bash batch_check.sh %s %s" % (program, " ".join("%d %d %d stars_check_%s" % (s + (name,)) for s, name in jobs)))
    print("      [batch_check.sh: %.0f s%s]" % (time.perf_counter() - t1, "; it said: " + out.strip()[:300] if out.strip() else ""), flush=True)
    bad, again_bad, count, threads = 0, 0, 0, (1, 2, 3, 4, 8)
    for s, name in jobs:
        stem = os.path.join(cdir, "batch_check_out", "%s_%d_%d_%d_stars_check_%s" % ((program,) + s + (name,)))
        n = len(sets[name])
        try:
            alone = read_lines(stem + "_alone.txt")
            with open(stem + "_alone.bin", "rb") as fh:
                alone_bin = fh.read()
            runs = {}
            for tag in ["b%d" % t for t in threads] + ["s%d" % t for t in threads] + ["again"]:
                with open(stem + "_%s.bin" % tag, "rb") as fh:
                    runs[tag] = (read_lines(stem + "_%s.txt" % tag), fh.read())
        except (OSError, ValueError) as exc:
            print("      %s at %s, %s: no output (%s)" % (program, s, name, exc))
            bad += 3 * n
            continue
        size = len(runs["b1"][1]) // n if n else 0
        expected, k0 = [], 0                         # for each star: (value bits, code, dump bytes), as it came out alone
        for line in alone:
            if line[0] == "refused":
                expected.append((-1.0e300, line[1], bytes(size)))
            else:
                expected.append((float.fromhex(line[0]), "0", alone_bin[k0 * size:(k0 + 1) * size]))
                k0 += 1
        print("      %s at %s, %d stars (%s): refused %d, early exit (ln L = -inf) %d" % (program, s, n, name,
              sum(1 for e in expected if e[1] != "0"), sum(1 for e in expected if e[1] == "0" and e[0] < -1.0e299)))
        for tag, (lines, dump) in runs.items():
            for i in range(n):
                j = int(np.where(order[name] == i)[0][0]) if tag.startswith("s") else i      # where star i stands in this run
                v, c, d = expected[i]
                got = lines[j] if j < len(lines) else ["none", "none"]
                w = -1.0e300 if got[0] == "-inf" else (float(got[0]) if got[0] != "none" else math.nan)
                same = len(alone) == n and len(lines) == n and math.copysign(1, v) == math.copysign(1, w) and v.hex() == w.hex() \
                    and got[1] == c and dump[j * size:(j + 1) * size] == d
                count += 1
                if not same:
                    bad += 1
                    print("      %s at %s, %s, run %s: star %d differs from the same star alone" % (program, s, name, tag, i))
        if runs["again"] != runs["b8"]:
            again_bad += 1
    check("%s: batch against alone, every bit of ln L, the code and every bin of %d star-runs (1, 2, 3, 4, 8 threads; also shuffled), stars that differ"
          % (program, count), bad, 0)
    check("%s: the same batch twice on 8 threads, sets whose bits differ" % program, again_bad, 0)

# 7: the shared library through Python's ctypes, against the batch mode's own lines
if "nicer_fast" in wanted:
    print("\n--- libnicer_fast.so through nicer_fast_lib.py", flush=True)
    stem = os.path.join(cdir, "batch_check_out", "nicer_fast_10_16_24_stars_check_%s_b8.txt")
    differ = 0
    for name in ("sixty", "refused"):
        code, out = linux("python3 -I -B nicer_fast_lib.py tables.bin 10 16 24 8 stars_check_%s.txt" % name)
        try:
            with open(stem % name, encoding="ascii") as fh:
                want = fh.read().splitlines()
        except OSError:
            want = []
        got = out.splitlines()
        differ += sum(1 for i in range(max(len(want), len(got))) if i >= len(want) or i >= len(got) or want[i] != got[i])
        if code or not want:
            differ += 1
            print("      nicer_fast_lib.py on %s stars: exit code %d, %s" % (name, code, out.strip()[:300]))
    check("nicer_fast_lib.py (ctypes) on the sixty and the refused stars, 8 threads: lines that differ from the batch mode's", differ, 0)

    # 8: the race detector
    print("\n--- make race (the compiler's race detector, -fsanitize=thread)", flush=True)
    t1 = time.perf_counter()
    code, out = linux("make race", timeout=900)
    for line in out.splitlines()[-6:]:
        print("      " + line)
    print("      [make race: %.0f s]" % (time.perf_counter() - t1), flush=True)
    check("make race: 256 stars on 8 threads under the race detector, races reported (or the run failed)",
          out.count("WARNING: ThreadSanitizer") + (1 if code else 0), 0)

for path in [os.path.join(cdir, "stars_check_%s%s.txt" % (name, suffix)) for name in sets for suffix in ("", "_shuffled")]:
    if os.path.exists(path):
        os.remove(path)
shutil.rmtree(os.path.join(cdir, "batch_check_out"), ignore_errors=True)
print("\n%d checks: %d hold, %d fail%s ; %.0f s" % (len(results), sum(results), len(results) - sum(results),
      " ; not built or not run: %s" % ", ".join(missing) if missing else "", time.perf_counter() - T0))
sys.exit(2 if missing else (0 if all(results) else 1))
