# -*- coding: utf-8 -*-
"""Check of bend.py against the test article's own tables (Bogdanov et al. 2019, supplementary files):
bending angle and lensing factor, and travel-time delay, for the test star (1.4 Suns, 12 km).
usage: python -I -B -X utf8 test_bend.py        (the apjlab5968 folder is found through settings.cfg at the repository's top)"""
import math
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from bend import theta_tau, lensing_factor
import settings

ref = settings.ref_dir()
GM_SUN = 1.32712440018e20      # m^3/s^2 (the Sun's GM, an observed value)
C = 299792458.0
M, R = 1.4, 12.0e3
L = 2.0 * GM_SUN * M / C**2
u = L / R
print("L = %.3f m, u = L/R = %.6f" % (L, u))

d = np.loadtxt(os.path.join(ref, "deflection_angle_lensing_factor_test_CU.txt"))
psi_ref, alpha, lens_ref = d[:, 0], d[:, 1], d[:, 2]
x = np.cos(alpha)
th, ta = theta_tau(u, x)
lf = lensing_factor(u, x)
ok = alpha < 1.5607
print("bending angle: %d rows, alpha up to %.4f rad; largest relative difference %.2e (rows with alpha below 1.5607: %.2e)"
      % (len(d), alpha.max(), np.max(np.abs(th / psi_ref - 1.0)), np.max(np.abs(th[ok] / psi_ref[ok] - 1.0))))
print("lensing factor: largest relative difference %.2e (rows with alpha below 1.5607: %.2e)"
      % (np.max(np.abs(lf / lens_ref - 1.0)), np.max(np.abs(lf[ok] / lens_ref[ok] - 1.0))))

d = np.loadtxt(os.path.join(ref, "travel_time_delay_test_CU.txt"))
psi_ref, alpha, dt_ref = d[:, 0], d[:, 1], d[:, 2]
m = alpha > 0
th, ta = theta_tau(u, np.cos(alpha[m]))
dt = ta * R / C
print("delay: %d rows; largest relative difference %.2e; largest delay %.4e s (reference %.4e s)"
      % (m.sum(), np.max(np.abs(dt / dt_ref[m] - 1.0)), dt.max(), dt_ref.max()))

# self-consistency: more nodes, and Codex's relation c dT/dTheta = b
xs = np.array([1e-6, 1e-3, 0.05, 0.2, 0.25, 0.3, 0.7, 0.999])
a1 = theta_tau(u, xs, n=160); a2 = theta_tau(u, xs, n=640)
print("160 nodes against 640 nodes: Theta %.1e, tau %.1e" % (np.max(np.abs(a1[0] - a2[0])), np.max(np.abs(a1[1] - a2[1]))))
xg = np.linspace(0.02, 0.98, 49); hh = 1e-5
tp = theta_tau(u, xg + hh); tm = theta_tau(u, xg - hh)
kb = np.sqrt((1 - xg**2) / (1 - u))
print("c dT/dTheta against b (both in units of R): largest difference %.1e" % np.max(np.abs((tp[1] - tm[1]) / (tp[0] - tm[0]) - kb)))
