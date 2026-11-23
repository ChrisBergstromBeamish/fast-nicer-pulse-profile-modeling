# -*- coding: utf-8 -*-
"""Check of public_likelihood.py against the two numbers the public code ships for this step alone
(X-PSI, xpsi/tests/test_background_marginalisation.py: one channel, exposure 1.0e6 s, phase shift 0.51522312, 32 counts,
a 64-point signal; expected -4521.763801536406 with support [0, 1e-11] and -2772.718471175409 with no upper limit, both to
a relative tolerance of 1e-8). The inputs below are copied from that test file. X-PSI itself is not run."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import public_likelihood as pl

exposure_time = 1.0e6
phase_shift = 0.51522312
counts = [0, 0, 6, 131, 445, 1121, 1965, 2967, 3861, 4725, 5207, 5702, 5912, 6022,
          5651, 5293, 4959, 4516, 3855, 3407, 2900, 2377, 2005, 1601, 1195, 902, 599,
          334, 151, 58, 8, 0]
signal = [0.1312685255677862, 0.1244573842611507, 0.11735800771358744,
          0.110104659697357, 0.10280996148537738, 0.09556550157583055,
          0.0884432158487088, 0.08149750256166845, 0.07476764890001181,
          0.06827980979061828, 0.06204913824314168, 0.05608198755342002,
          0.05037773907206468, 0.04493038960172567, 0.03972988623083532,
          0.03476313858800521, 0.030014885102893008, 0.025468159584118875,
          0.021104738892653817, 0.016903992233303595, 0.012850280618456302,
          0.009153761613574292, 0.0060369913011612125, 0.003588428269002173,
          0.0018201316624806847, 0.0006985786687345404, 0.00013955225824213515,
          2.942156216152547e-06, 1.3666828401156535e-36, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 3.0117271162432395e-10, 5.655308922833262e-05,
          0.0009806827413526417, 0.0038040432854325266, 0.009142175813641658,
          0.017318757504263326, 0.028368022060331424, 0.041997897093318394,
          0.05742927955365583, 0.07312413392268027, 0.08819357752072775,
          0.10231526616772711, 0.11518864065189004, 0.12655365326319978,
          0.13620586164399945, 0.14400508855527985, 0.14988005012774142,
          0.15382729387100896, 0.15590563078819905, 0.15622702770198071,
          0.15494509787599256, 0.15224282805060138, 0.14832037071555773,
          0.14338387244889514, 0.13763594437424936, 0.1312685255677862]

d = np.array(counts, dtype=float)
rate = np.array([signal], dtype=float)
star = pl.phase_integrate(rate, np.linspace(0.0, 1.0, 64), np.linspace(0.0, 1.0, 33), phase_shift)[0]
for support, expected in (([0.0, 1.0e-11], -4521.763801536406), ([0.0, -1.0], -2772.718471175409)):
    lower, upper, B_for, B = pl.channel_limits(d, star, exposure_time, support)
    val = pl.marginal_loglike_channel(d, star, exposure_time, support)
    print("support %s: ours %.12f ; the public code's shipped value %.12f ; difference %.3e (relative %.2e) ; limits %.6e to %.6e"
          % (support, val, expected, val - expected, abs(val - expected) / abs(expected), lower, upper), flush=True)
