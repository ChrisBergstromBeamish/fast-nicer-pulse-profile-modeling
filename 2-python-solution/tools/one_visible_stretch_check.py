"""Is the visible part of a ring's turn always ONE stretch that starts at the point nearest the observer?
The fast kernel assumes it (one pair of end directions per ring). Codex's round-3 review asked for the check.
For many stars inside the public example's prior (and faster-spinning ones), and for many rings each, mu is worked out at
4000 points along the half-turn; the test counts how often visibility switches on again after it has switched off, and how
much of the ring's light such a second stretch could hold (the integral of mu over it against the integral over all)."""
import math, os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
from consts import C_CM_S, GM_SUN
rng = np.random.default_rng(21)
phi = np.linspace(0.0, math.pi, 4001)
n_rings_total, n_second, worst_share, worst_case = 0, 0, 0.0, None
xt = None
for case in range(400):
    M = rng.uniform(1.0, 2.2); Rk = rng.uniform(9.5, 15.0); nu = rng.choice([100.0, 314.0, 500.0, 716.0])
    incl = math.acos(rng.uniform(-1.0, 1.0))
    GM = GM_SUN * M; L = 2 * GM / C_CM_S**2; Req = Rk * 1e5; Om = 2 * math.pi * nu
    o2 = (Om**2 * Req**3 / GM) * (-0.788 + 1.030 * GM / (Req * C_CM_S**2))
    if not (0.13 < L / Req and L / (Req * (1 + o2)) < 0.56):
        continue
    theta = np.arccos(rng.uniform(-1.0, 1.0, 40))
    ct, st = np.cos(theta), np.sin(theta)
    R = Req * (1 + o2 * ct * ct); u = L / R; f = np.sqrt(1 - u)
    h = -2 * Req * o2 * ct * st / (R * f); rt = np.sqrt(1 + h * h); cg, sg = 1 / rt, h / rt
    ci, si = math.cos(incl), math.sin(incl)
    for k in range(theta.size):
        x_floor = -(abs(sg[k]) + 0.02)
        xs = x_floor + (1.0 - x_floor) * np.linspace(0.0, 1.0, 1500) ** 1.5
        xs[-1] = 1.0 - 1e-13
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            th = fk.bend(np.full(xs.size, u[k]), xs, n=32, need_delay=False)[0]
        if not np.all(np.isfinite(th)) or th[0] >= math.pi:
            continue
        C = ci * ct[k] + si * st[k] * np.cos(phi)
        psi = np.arccos(np.clip(C, -1.0, 1.0))
        reach = psi < th[0]
        x = np.interp(psi, th[::-1], xs[::-1])
        D = (ci - ct[k] * C) / (st[k] * np.sqrt(np.maximum(1.0 - C * C, 1e-300)))
        mu = x * cg[k] + np.sqrt(np.maximum(1.0 - x * x, 0.0)) * sg[k] * np.clip(D, -1.0, 1.0)
        vis = reach & (mu > 0.0)
        n_rings_total += 1
        sw = np.diff(vis.astype(int))
        ons = np.where(sw == 1)[0]
        if vis[0] and ons.size == 0:
            continue
        if not vis.any():
            continue
        # a stretch that does not start at phi = 0
        first_off = np.where(~vis)[0][0] if (~vis).any() else vis.size
        second = vis.copy(); second[:first_off] = False
        share = float(np.sum(np.where(second, mu, 0.0)) / max(np.sum(np.where(vis, mu, 0.0)), 1e-300))
        if second.any():
            n_second += 1
            if share > worst_share:
                worst_share, worst_case = share, (M, Rk, nu, math.degrees(incl), math.degrees(theta[k]), float(sg[k]), int(second.sum()), float(mu[second].max()))
print("%d rings on %d stars: rings with a visible stretch that does not start at the nearest point: %d" % (n_rings_total, 400, n_second))
if worst_case:
    print("largest share of a ring's mu-weighted turn in such a stretch: %.2e  (star %.2f Suns, %.2f km, %.0f Hz, observer at %.1f deg, ring at %.1f deg, sin(tilt) %.4f; %d of 4001 points, largest mu there %.2e)"
          % ((worst_share,) + worst_case))
