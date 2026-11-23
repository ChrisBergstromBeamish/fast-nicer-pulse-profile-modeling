# -*- coding: utf-8 -*-
"""One command that repeats the checks of the Python solution and says for each whether it holds.

usage, from the repository's top:   python -I -B -X utf8 2-python-solution/check.py [--slow]
  without --slow   the quick checks (about ten seconds on one processor core of a laptop)
  with --slow      also the slow checks (about ten minutes): the approximation of the path integrals against the integrals,
                   the sixteen agreed reference curves, the 247 hard made-up channels, the public code's second check point,
                   the plain way on 36 round stars, the one-stretch test, the shipped counts file, counting at emission against
                   the usual way, and the plain way at fine patches on the quick example's star and the second check star.
                   (Two slower comparisons keep their own commands: tools/kernel_against_plain_star_many.py, about fifteen
                   minutes, and tools/kernel_against_plain_star.py, about six minutes.)
Each line says: ok or FAIL, what is compared, the difference found, the limit. The last line counts them. The exit code is 0 if
every check holds, 1 if any fails, 2 if the settings or the reference material are missing.
The compiled solution is not needed here; it has its own command, 3-compiled-solution/check.py.
The reference material is found through settings.cfg at the repository's top (1-setup-and-documentation/README.md).
Every number that this command checks has a row in 1-setup-and-documentation/CLAIMS.md.
"""
import math
import os
import re
import subprocess
import sys
import time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "code"))
sys.path.insert(0, os.path.join(HERE, "tools"))
import numpy as np
import settings

try:
    data_dir, ref_dir = settings.data_dir(), settings.ref_dir()
except settings.SettingsError as exc:
    print("cannot run the checks: %s" % exc)
    sys.exit(2)
slow = "--slow" in sys.argv[1:]
unknown = [a for a in sys.argv[1:] if a != "--slow"]
if unknown:
    print("usage: python -I -B -X utf8 2-python-solution/check.py [--slow]")
    sys.exit(2)
import fast_kernel as fk                    # noqa: E402  (after the settings, so that a missing settings file is reported first)
import plain_star as pst                    # noqa: E402
import public_likelihood as pl              # noqa: E402
import quick_example as qe                  # noqa: E402
from bend import theta_tau                  # noqa: E402
from consts import BB_CONST, PC_CM          # noqa: E402

results = []
T0 = time.perf_counter()
print("Checks of the Python solution%s. The public example's data: %s ; the agreed curves and tables: %s"
      % (" (quick and slow)" if slow else " (quick)", data_dir, ref_dir), flush=True)


def check(name, value, limit):
    """One line: does the difference found stay within the limit?"""
    ok = bool(value <= limit)
    results.append(ok)
    print("%s  %-104s %.2e (limit %.3g)" % ("ok  " if ok else "FAIL", name, value, limit), flush=True)


# ---- the quick checks ---------------------------------------------------------------------------------------------------
# 1. light's path against the test article's tables (Bogdanov et al. 2019, supplementary files; the star of 1.4 Suns and 12 km)
cu = np.loadtxt(os.path.join(ref_dir, "deflection_angle_lensing_factor_test_CU.txt"))
L_over_R = 2.0 * 1.4 * 1.32712440018e26 / (2.99792458e10 ** 2) / 12.0e5
alpha = cu[:, 1]                       # the table's columns: swept angle, leaving angle, lensing factor
sel = alpha < 1.5607
th = fk.bend(np.full(sel.sum(), L_over_R), np.cos(alpha[sel]), n=96, need_delay=False)[0]
check("swept angle against the article's table (its printed precision)", float(np.max(np.abs(th / cu[sel, 0] - 1.0))), 2e-5)
dl = np.loadtxt(os.path.join(ref_dir, "travel_time_delay_test_CU.txt"))
m = dl[:, 1] > 0.0
tau = theta_tau(L_over_R, np.cos(dl[m, 1]))[1] * 12.0e5 / 2.99792458e10
check("delay behind a radial ray against the article's table (its printed precision)", float(np.max(np.abs(tau / dl[m, 2] - 1.0))), 2e-5)
xg, hh = np.linspace(0.02, 0.98, 49), 1e-5
tp, tm = theta_tau(L_over_R, xg + hh), theta_tau(L_over_R, xg - hh)
check("the delay's rate of change with the swept angle against b (an identity of the two integrals)",
      float(np.max(np.abs((tp[1] - tm[1]) / (tp[0] - tm[0]) - np.sqrt((1.0 - xg ** 2) / (1.0 - L_over_R))))), 1e-7)
# 2. the likelihood step against the public code's two shipped unit-test values (xpsi/tests/test_background_marginalisation.py)
exposure_time, phase_shift = 1.0e6, 0.51522312
cnt = np.array([0, 0, 6, 131, 445, 1121, 1965, 2967, 3861, 4725, 5207, 5702, 5912, 6022, 5651, 5293, 4959, 4516, 3855, 3407, 2900, 2377, 2005,
                1601, 1195, 902, 599, 334, 151, 58, 8, 0], dtype=float)
sig = np.array([0.1312685255677862, 0.1244573842611507, 0.11735800771358744, 0.110104659697357, 0.10280996148537738, 0.09556550157583055,
                0.0884432158487088, 0.08149750256166845, 0.07476764890001181, 0.06827980979061828, 0.06204913824314168, 0.05608198755342002,
                0.05037773907206468, 0.04493038960172567, 0.03972988623083532, 0.03476313858800521, 0.030014885102893008, 0.025468159584118875,
                0.021104738892653817, 0.016903992233303595, 0.012850280618456302, 0.009153761613574292, 0.0060369913011612125,
                0.003588428269002173, 0.0018201316624806847, 0.0006985786687345404, 0.00013955225824213515, 2.942156216152547e-06,
                1.3666828401156535e-36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0117271162432395e-10, 5.655308922833262e-05,
                0.0009806827413526417, 0.0038040432854325266, 0.009142175813641658, 0.017318757504263326, 0.028368022060331424,
                0.041997897093318394, 0.05742927955365583, 0.07312413392268027, 0.08819357752072775, 0.10231526616772711, 0.11518864065189004,
                0.12655365326319978, 0.13620586164399945, 0.14400508855527985, 0.14988005012774142, 0.15382729387100896, 0.15590563078819905,
                0.15622702770198071, 0.15494509787599256, 0.15224282805060138, 0.14832037071555773, 0.14338387244889514, 0.13763594437424936,
                0.1312685255677862])
st = pl.phase_integrate(sig[None, :], np.linspace(0.0, 1.0, 64), np.linspace(0.0, 1.0, 33), phase_shift)[0]
for sup, expected in (([0.0, 1.0e-11], -4521.763801536406), ([0.0, -1.0], -2772.718471175409)):
    v = pl.marginal_loglike_channel(cnt, st, exposure_time, sup)
    check("likelihood step against the public code's shipped unit value %.6f (relative)" % expected, abs(v / expected - 1.0), 1e-8)
# 3. the quick example
counts, bg, support = qe.load_data(data_dir)
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph = np.linspace(0.0, 1.0, qe.NB + 1)


def like(P, sup, nr, n, nE, star_only=False):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    M, D, kT, inc = qe.star_args(P)
    S = fk.spot_counts(M, P["radius"], P["freq"], D, inc, P["colat"], P["ang_radius"], kT, E, ph, phase0=0.5 - P["phase_shift"], n_rings=nr, n=n)
    star = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * (W @ S.T)
    return star if star_only else pl.marginal_loglike_fast(counts, star, qe.T_EXP, sup, fact)


v_conv = like(qe.P, support, 48, 32, 48)
v_work = like(qe.P, support, 10, 16, 24)
check("quick example: ours (converged) against the public code's shipped likelihood; its own tolerance", abs(v_conv - qe.EXPECTED), qe.TOL)
check("quick example: ours at the working setting against ours converged", abs(v_work - v_conv), 0.02)
star = like(qe.P, support, 10, 16, 24, star_only=True)
check("quick example: the fast likelihood form against the step-by-step form", abs(pl.marginal_loglike_fast(counts, star, qe.T_EXP, support, fact)
                                                                                 - pl.marginal_loglike(counts, star, qe.T_EXP, support)), 1e-8)
# 4. the fast kernel against the plain way on the example's star (coarse patches, four energies)
E4 = np.array([0.2, 0.5, 1.0, 1.6])
M, D, kT, inc = qe.star_args(qe.P)
NB, SUB = 32, 16
mid = (np.arange(NB * SUB) + 0.5) / (NB * SUB)
Sk = fk.spot_counts(M, 12.0, 314.0, D, inc, qe.P["colat"], qe.P["ang_radius"], kT, E4, np.linspace(0.0, 1.0, NB + 1), n_rings=24, n=24)
Sp = pst.light_curve(M, 12.0, 314.0, D, inc, qe.P["colat"], qe.P["ang_radius"], kT, E4, mid, n_zeta=32, n_eta=128, n_phi=1024).reshape(NB, SUB, 4).mean(axis=1) / NB
check("fast kernel against the plain way, bin by bin, of the mean bin (plain way at coarse patches)", float(np.max(np.abs(Sk - Sp)) / Sp.mean()), 2e-4)
# 5. one agreed reference curve, directly
d = np.loadtxt(os.path.join(ref_dir, "OS1e_test_IM.txt"))
w = 1.0 / 2048.0
edges = np.empty(2 * d.shape[0])
edges[0::2], edges[1::2] = d[:, 0] - 0.5 * w, d[:, 0] + 0.5 * w
y = fk.spot_counts(1.4, 12.0, 600.0, 200.0 * PC_CM, math.radians(30.0), math.radians(60.0), 1.0, 0.35, [1.0], edges, n_rings=24, n=24)[0::2, 0] / w
check("fast kernel against the agreed reference curve OS1e, of its median (the article's target 1e-3)", float(np.max(np.abs(y - d[:, 1])) / np.median(d[:, 1])), 1e-3)
n_quick = len(results)

# ---- the slow checks: each is its own script in tools/, run as a command; its last line "RESULT key=value ..." is read here --
SLOW = (
    ("path_approximation_check.py", [], (
        ("swept_angle", "the approximation of the path integrals against the integrals themselves on 6000 random directions: swept angle, rad", 1e-11),
        ("delay", "the approximation of the path integrals against the integrals themselves on 6000 random directions: delay, in units of R/c", 1e-11))),
    ("kernel_against_reference_curves.py", [], (
        ("extended_24_24", "the twelve agreed reference curves with extended spots, 24 rings per panel and 24 points, of the median (the article's target)", 1e-3),
        ("extended_10_16", "the same twelve curves at the working setting (10 rings per panel, 16 points), of the median", 1e-3),
        ("pointlike_24_24", "the four curves with point-like spots (radius 0.01), of the median: above the article's target at eclipse edges", 1e-2))),
    ("fast_likelihood_stress.py", [], (
        ("largest", "the fast likelihood form against the step-by-step form on 247 hard made-up channels, per channel", 1e-8),)),
    ("second_check_point.py", [], (
        ("converged_minus_shipped_vs_found_before", "the public code's second check: ours converged minus its shipped value, against the +6.518 found before (ours FAILS that check)", 5e-3),)),
    ("kernel_against_plain_way.py", [], (
        ("flagged", "round stars, 36 geometries against the plain way: geometries above 2e-4 of the largest bin (two known, with a sliver of the spot in view)", 2),
        ("largest", "round stars, 36 geometries against the plain way: largest difference of the largest bin", 4e-4))),
    ("one_visible_stretch_check.py", [], (
        ("second_stretch_rings", "rings whose visible stretch does not start at the nearest point, of 14968 rings on 400 stars (13 found before; the kernel takes one stretch)", 13),)),
    ("quick_example_counts_check.py", [], (
        ("total_ratio_abs", "quick example: our total expected counts against the shipped counts file, |(ours / file) - 1| (1.45e-5 found before)", 5e-5),
        ("interior_channels_largest", "quick example: per channel summed over the turn, channels 2 to 289, against the shipped counts file (1.8e-5 found before)", 5e-5),
        ("largest_rel", "quick example: largest single difference against the shipped counts file (1.16e-3 found before, in the edge channels and faint bins)", 2e-3))),
    ("arrival_grid_check.py", [], (
        ("lnL_129", "the even grid that locates arrival times: 129 points against 2049, the quick example's ln L at the working setting", 1e-5),
        ("bins_129", "the even grid that locates arrival times: 129 points against 2049, the bins, of the mean bin", 1e-8))),
    ("count_at_emission_check.py", [], (
        ("largest_64", "counting photons at emission against the usual way (400,000 phases), Fourier coefficients of one element's light curve, 64 nodes", 1e-8),)),
    ("likelihood_from_plain_star.py", [], (
        ("quick_96x384_bins", "quick example's star: fast kernel (48 rings per panel, 32 points) against the plain way at 96 x 384 patches, bin by bin, of the mean bin", 1e-5),
        ("quick_96x384_minus_kernel_abs", "quick example's star: the plain way's ln L (96 x 384 patches) minus the fast kernel's, +0.0009 found before", 5e-3),
        ("second_96x384_bins", "second check star: fast kernel against the plain way at 96 x 384 patches, bin by bin, of the mean bin", 1e-5),
        ("second_96x384_minus_kernel_abs", "second check star: the plain way's ln L (96 x 384 patches) minus the fast kernel's, -0.089 found before (both above the shipped value)", 0.2))),
)
# counts that the wording of the slow checks relies on: each becomes a check that the script worked on exactly that many
COUNTS = {"path_approximation_check.py": (("directions", 6000, "random directions in the approximation check"),),
          "fast_likelihood_stress.py": (("channels", 247, "hard made-up channels kept (of 271 made; the rest the public code's early exit)"),),
          "kernel_against_plain_way.py": (("geometries", 36, "round-star geometries"), ("refused", 0, "round-star geometries refused")),
          "one_visible_stretch_check.py": (("rings", 14968, "rings in the one-stretch test"), ("stars", 400, "stars in the one-stretch test"))}
if slow:
    print("\nthe slow checks (each is a script in 2-python-solution/tools, run here as a command; its own output is shown indented):", flush=True)
    for script, args, items in SLOW:
        t1 = time.perf_counter()
        cmd = [sys.executable, "-I", "-B", "-X", "utf8", os.path.join(HERE, "tools", script)] + args
        done = subprocess.run(cmd, capture_output=True, text=True)
        out = (done.stdout or "") + (done.stderr or "")
        for line in out.splitlines():
            print("      " + line)
        print("      [%s: %.0f s]" % (script, time.perf_counter() - t1), flush=True)
        found = {}
        for line in out.splitlines():
            if line.startswith("RESULT "):
                for key, value in re.findall(r"(\w+)=(\S+)", line):
                    found[key] = float(value)
        for key, name, limit in items:
            if key in found:
                check(name, found[key], limit)
            else:
                check(name + " [the script gave no %s]" % key, math.inf, limit)
        for key, expected, what in COUNTS.get(script, ()):
            check("%s: the number of %s is %d, as the wording above says" % (script, what, expected), abs(found.get(key, math.inf) - expected), 0)
print("\n%d quick checks: %d hold, %d fail%s ; %.0f s" % (n_quick, sum(results[:n_quick]), n_quick - sum(results[:n_quick]),
      "" if not slow else " ; %d slow checks: %d hold, %d fail" % (len(results) - n_quick, sum(results[n_quick:]), len(results) - n_quick - sum(results[n_quick:])),
      time.perf_counter() - T0))
sys.exit(0 if all(results) else 1)
