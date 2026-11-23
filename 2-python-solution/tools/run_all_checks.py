# -*- coding: utf-8 -*-
"""One command that repeats the quick checks of the NICER programs and says for each whether it still holds.
(The slow checks are not in here: the plain way on many geometries, the sixty prior draws, the timing. See the other files
in this folder.)

usage: python -I -B -X utf8 run_all_checks.py <the public example's Data folder> <folder of the agreed reference curves>
"""
import math
import os
import subprocess
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
import plain_star as pst
import public_likelihood as pl
import quick_example as qe
from consts import BB_CONST, PC_CM

data_dir, ref_dir = sys.argv[1], sys.argv[2]
results = []


def check(name, value, limit, unit=""):
    ok = value <= limit
    results.append(ok)
    print("%s  %-92s %.2e (limit %.0e)%s" % ("ok  " if ok else "FAIL", name, value, limit, unit), flush=True)


# 1. light's path against the test article's tables
cu = np.loadtxt(os.path.join(ref_dir, "deflection_angle_lensing_factor_test_CU.txt"))
L_over_R = 2.0 * 1.4 * 1.32712440018e26 / (2.99792458e10 ** 2) / 12.0e5
alpha = cu[:, 1]                       # the table's columns: swept angle, leaving angle, lensing factor
sel = alpha < 1.5607
th = fk.bend(np.full(sel.sum(), L_over_R), np.cos(alpha[sel]), n=96, need_delay=False)[0]
check("swept angle against the article's table (its printed precision)", float(np.max(np.abs(th / cu[sel, 0] - 1.0))), 2e-5)
# 2. the likelihood step against the public code's two shipped unit-test values
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
# 6. the compiled copy, if it is built
wsl_c = (lambda p: "/mnt/" + p[0].lower() + p[2:].replace(chr(92), "/") if p[1:2] == ":" else p)(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")))  # the same folder as the Linux side of a Windows machine sees it
cmd = "cd %s && ./nicer_fast tables.bin 10 16 24 0 1.4 12 1.0 %r 0.0 %r 0.75 6.7" % (wsl_c, qe.P["cos_i"], qe.P["colat"])
try:
    out = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", cmd], capture_output=True, text=True, timeout=120).stdout.replace(chr(0), "")
    line = [l for l in out.splitlines() if l.startswith("ln L =")]
    if line:
        check("compiled copy against the Python at the working setting (ln L)", abs(float(line[0].split("=")[1]) - v_work), 1e-6)
    else:
        print("skip  the compiled copy gave no answer (not built?)")
except Exception as exc:                                               # noqa: BLE001
    print("skip  the compiled copy could not be run (%s)" % type(exc).__name__)
print("\n%d checks, %d hold, %d fail" % (len(results), sum(results), len(results) - sum(results)))
