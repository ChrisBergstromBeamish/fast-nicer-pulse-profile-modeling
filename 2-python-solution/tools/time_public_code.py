# -*- coding: utf-8 -*-
"""Stopwatch for the public code X-PSI on its own shipped quick example. Runs INSIDE X-PSI's own environment, once it is
installed (written on 22 Nov 2026, before the install was allowed; first run on 23 Nov 2026).

It builds the objects exactly as examples/examples_fast/Modules/main.py does (the lines are copied from it), at the
numerical settings given on the command line, then
  1. evaluates the likelihood at the example's star p and prints it beside the shipped value -32552.616905;
  2. times the likelihood at points near p, a new point at every call (as a sampler does), caches cleared.
One line of JSON is appended to the output file.

Known limit, found when it was first run (23 Nov 2026), for the second check: the points near the star are the star with
each number moved by up to 1 part in 1000 either way, and that star's radius, 10, is the public code's lower limit of the
radius. So about half of the points lie below that limit and are refused at once (the public code answers with no
evaluation); their times are not times of an evaluation. compare_with_public_code.py leaves the refused points out and
gives times and ratios for the quick check only. The value and the bins at the star itself are not affected.

usage: python time_public_code.py <xpsi source root> <output .jsonl> [sqrt_num_cells num_leaves num_rays num_energies [n_calls [quick|second]]]
  second: the public code's second shipped check instead (xpsi/tests/test_likelihood.py: radius 10, log10 T 6.8, no limits
  on the background; shipped value -47881.27817666349 at 128 energies).
"""
import json
import math
import os
import sys
import time

import numpy as np

root = sys.argv[1]
out_path = sys.argv[2]
sqrt_num_cells = int(sys.argv[3]) if len(sys.argv) > 3 else 32
num_leaves = int(sys.argv[4]) if len(sys.argv) > 4 else 64
num_rays = int(sys.argv[5]) if len(sys.argv) > 5 else 512
num_energies = int(sys.argv[6]) if len(sys.argv) > 6 else 64
n_calls = int(sys.argv[7]) if len(sys.argv) > 7 else 12
mode = sys.argv[8] if len(sys.argv) > 8 else "quick"
quick = (sqrt_num_cells, num_leaves, num_rays, num_energies) == (32, 64, 512, 64)

import xpsi

modules = os.path.join(root, "examples", "examples_fast", "Modules")
sys.path.append(modules)
from CustomInstrument import CustomInstrument
from CustomSignal import CustomSignal
from CustomPhotosphere import CustomPhotosphere
from CustomPrior import CustomPrior

# ---- copied from main.py (data, instrument, support) ----
data_loaded = np.loadtxt(os.path.join(modules, "..", "Data", "xpsi_good_realisation.dat"), dtype=np.double)
data = xpsi.Data(data_loaded,
                 channels=np.arange(10, 301),
                 phases=np.linspace(0.0, 1.0, 33),
                 first=0,
                 last=290,
                 exposure_time=1000.0)
channel_number = np.arange(0, 1501)
energy_low = np.arange(0, 15.01, 0.01)
energy_high = energy_low + 0.01
channel_edges = np.array([list(channel_number), list(energy_low), list(energy_high)]).T
arf_energy_low = [0.1]
arf_energy_high = [0.105]
arf_val = [1800]
counter = 1
while arf_energy_low[-1] <= 14.995:
    arf_energy_low.append(arf_energy_low[-1] + 0.005)
    arf_energy_high.append(arf_energy_high[-1] + 0.005)
    arf_val.append(1800)
    counter += 1
ARF = np.array([list(arf_energy_low), list(arf_energy_high), list(arf_val)]).T
RMF = np.diag(np.full(counter, 1))
Instrument = CustomInstrument.from_response_files(ARF=ARF, RMF=RMF, channel_edges=channel_edges,
                                                  max_input=301, min_input=10, channel=[10, 301])
data_background = np.loadtxt(os.path.join(modules, "..", "Data", "xpsi_synthetic_background.dat"))
registered_true_background_rate = data_background[:, 1] / data.exposure_time
support = np.ones((len(registered_true_background_rate), 2), dtype=np.double)
support[:, 0] = 0.9 * registered_true_background_rate
support[:, 1] = 1.1 * registered_true_background_rate

if mode == "second":           # as in xpsi/tests/test_likelihood.py: no limits on the background
    signal = CustomSignal(data=data, instrument=Instrument, interstellar=None, cache=True,
                          workspace_intervals=1000, epsrel=1.0e-8, epsilon=1.0e-3, sigmas=10.0)
else:
    signal = CustomSignal(data=data, instrument=Instrument, support=support, interstellar=None, cache=True,
                          workspace_intervals=1000, epsrel=1.0e-8, epsilon=1.0e-3, sigmas=10.0)
bounds = dict(distance=(0.5, 2), mass=(1.0, 1.6), radius=(10, 13), cos_inclination=(0, 1))
spacetime = xpsi.Spacetime(bounds, values=dict(frequency=314.0))
bounds = dict(super_colatitude=(0.001, math.pi / 2 - 0.001),
              super_radius=(0.001, math.pi / 2 - 0.001),
              phase_shift=(-0.25, 0.75),
              super_temperature=(6., 7.))
hot_spot = xpsi.HotRegion(bounds=bounds, values={}, symmetry=True, omit=False, cede=False, concentric=False,
                          sqrt_num_cells=sqrt_num_cells,
                          min_sqrt_num_cells=(16 if quick else sqrt_num_cells // 2),
                          max_sqrt_num_cells=(64 if quick else 4 * sqrt_num_cells),
                          num_leaves=num_leaves, num_rays=num_rays,
                          is_antiphased=True, image_order_limit=3, prefix='p')
photosphere = CustomPhotosphere(hot=hot_spot, elsewhere=None, values=dict(mode_frequency=spacetime['frequency']))
star = xpsi.Star(spacetime=spacetime, photospheres=photosphere, name='algendy')
prior = CustomPrior()
likelihood = xpsi.Likelihood(star=star, signals=signal, num_energies=num_energies, threads=1,
                             externally_updated=False, prior=prior)
# ---- end of the copied part ----

p = [1.4, 12, 1., math.cos(60 * np.pi / 180), 0.0, 70 * np.pi / 180, 0.75, 6.7]
EXPECTED = -3.2552616905e+04
if mode == "second":
    p = [1.4, 10, 1., math.cos(60 * np.pi / 180), 0.0, 70 * np.pi / 180, 0.75, 6.8]
    EXPECTED = -47881.27817666349
t0 = time.perf_counter()
v0 = float(likelihood(p, force=True))
t_first = time.perf_counter() - t0
print("%s check ; settings: sqrt_num_cells %d, leaves %d, rays %d, energies %d" % (mode, sqrt_num_cells, num_leaves, num_rays, num_energies))
print("at p: ln L = %.9f ; shipped value %.9f ; difference %+.6f ; first call %.4f s" % (v0, EXPECTED, v0 - EXPECTED, t_first), flush=True)

# What the public code worked out at p, saved for a bin-by-bin comparison with ours (each item only if the object has it).
saved = {}
for name, get in (("rate", lambda: np.asarray(signal.signals[0])),                      # counts/s, (channels, signal phases)
                  ("expected_counts", lambda: np.asarray(signal.expected_counts)),       # source plus best background, (channels, bins)
                  ("background_counts", lambda: np.asarray(signal.background_signal)),   # best background, counts over the exposure
                  ("energies", lambda: np.asarray(signal.energies)),
                  ("signal_phases", lambda: np.asarray(signal.phases[0])),
                  ("cell_area", lambda: np.asarray(hot_spot._super_cellArea)),
                  ("cell_colatitude", lambda: np.asarray(hot_spot._super_theta)),
                  ("cell_azimuth", lambda: np.asarray(hot_spot._super_phi))):
    try:
        saved[name] = get()
    except Exception as exc:                                                             # noqa: BLE001 (a missing attribute is not fatal)
        print("not saved: %s (%s)" % (name, type(exc).__name__))
npz_path = "%s.%s_cells%d_leaves%d_rays%d_energies%d.npz" % (out_path, mode, sqrt_num_cells, num_leaves, num_rays, num_energies)
np.savez_compressed(npz_path, value_at_p=v0, **saved)
print("saved: %s (%s)" % (npz_path, ", ".join("%s %s" % (k, v.shape) for k, v in saved.items())), flush=True)

rng = np.random.default_rng(20261122)
scale = np.array([1e-3, 1e-3, 1e-3, 1e-3, 0.0, 1e-3, 1e-3, 1e-4])
points, values, times = [], [], []
for k in range(n_calls):
    q = list(np.array(p) * (1.0 + scale * rng.uniform(-1.0, 1.0, 8)))
    q[4] = 1e-3 * rng.uniform(-1.0, 1.0)
    t0 = time.perf_counter()
    v = float(likelihood(q, force=True))
    dt = time.perf_counter() - t0
    points.append([float(x) for x in q])
    values.append(v)
    times.append(dt)
    print("call %2d: ln L = %.6f ; %.4f s" % (k, v, dt), flush=True)
times_sorted = sorted(times)
rec = dict(mode=mode, star=[float(x) for x in p], settings=dict(sqrt_num_cells=sqrt_num_cells, num_leaves=num_leaves, num_rays=num_rays, num_energies=num_energies),
           xpsi_version=getattr(xpsi, "__version__", "?"), value_at_p=v0, shipped=EXPECTED, first_call_s=t_first,
           points=points, values=values, times_s=times,
           fastest_s=times_sorted[0], median_s=times_sorted[len(times_sorted) // 2],
           numpy=np.__version__, python=sys.version.split()[0])
with open(out_path, "a", encoding="utf-8") as fh:
    fh.write(json.dumps(rec) + "\n")
print("fastest %.4f s ; median %.4f s per likelihood evaluation" % (rec["fastest_s"], rec["median_s"]))
print("TIMING_DONE")
