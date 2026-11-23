# First check of one planned optimization (22 Nov 2026): count the photons at emission.
# The usual way, for each rotation phase: find the emission angle that reaches the observer (an inversion), then multiply by the
# lensing factor d(cos alpha)/d(cos psi) (a derivative). The planned way: integrate over the emission angle itself. Then the
# lensing factor cancels against the change of variable, and no inversion is needed. Both ways are computed here for one small
# surface element and compared through the Fourier coefficients of its light curve over one rotation.
#
# Light's path used here comes from clocks and rulers only (units: the length L = 2GM/c^2 is 1): clock rate f(r)^2 = 1 - 1/r,
# radial ruler dl = dr/f, light at speed c in every frame at rest, and the constant b = r sin(alpha)/f(r) along the path
# (the angular momentum law, accepted by the author as a premise on 23 Nov 2026). That gives d(psi)/dr = b / (r^2 sqrt(1 - b^2 f^2 / r^2)).
# Slow rotation and an isotropic surface are taken, so only the light path is tested.
import math
import numpy as np

_G, _W = np.polynomial.legendre.leggauss(600)            # computed once
_V = 0.5 * (_G + 1.0); _WV = 0.5 * _W                    # v = sqrt(1 - R/r) in [0, 1]
_WW = 1.0 - _V * _V                                      # R/r

def psi_forward(u, x):
    """Angle psi swept by the ray, from cos(alpha) = x at the surface (radius 1/u) out to far away. Forward, no inversion."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    s2 = (1.0 - x * x) / (1.0 - u)                       # (b/R)^2
    core = (_WW * _WW * (1.0 - u * _WW))[None, :]
    out = np.empty(x.size)
    for a in range(0, x.size, 2000):                     # in pieces, to keep the memory small
        s = s2[a:a + 2000, None]
        den = np.maximum(1.0 - s * core, 1e-300)
        out[a:a + 2000] = (2.0 * np.sqrt(s) * _V[None, :] / np.sqrt(den)) @ _WV
    return out

def usual_way(u, inc, theta, mmax, nphi=400000, ntab=8001):
    xt = np.cos(np.linspace(0.0, math.pi / 2, ntab))     # table of cos(alpha), from 1 down to 0
    ct = np.cos(psi_forward(u, xt))                      # table of cos(psi)
    dxdc = np.gradient(xt, ct)                           # lensing factor, by differences of the table
    phi = (np.arange(nphi) + 0.5) * 2 * math.pi / nphi
    c = math.cos(inc) * math.cos(theta) + math.sin(inc) * math.sin(theta) * np.cos(phi)
    vis = c >= ct[-1]
    x = np.interp(c, ct[::-1], xt[::-1])                 # the inversion
    D = np.interp(c, ct[::-1], dxdc[::-1])
    g = np.where(vis, x * D, 0.0)
    return np.array([np.mean(g * np.cos(m * phi)) for m in range(mmax + 1)])

def at_emission(u, inc, theta, mmax, n):
    c1 = math.cos(inc - theta); c2 = math.cos(inc + theta)       # largest and smallest cos(psi) over one rotation
    def x_of_c(c):                                               # cos(alpha) at one end of the range: two numbers per ring
        if math.cos(float(psi_forward(u, [0.0])[0])) > c:
            return 0.0                                           # that end is hidden: the range stops at cos(alpha) = 0
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if math.cos(float(psi_forward(u, [mid])[0])) > c: hi = mid
            else: lo = mid
        return 0.5 * (lo + hi)
    x1 = x_of_c(c1); x2 = x_of_c(c2)
    t = (np.arange(n) + 0.5) * math.pi / n
    x = 0.5 * (x1 + x2) + 0.5 * (x1 - x2) * np.cos(t)            # nodes in cos(alpha)
    c = np.cos(psi_forward(u, x))                                # forward only
    cphi = np.clip((c - math.cos(inc) * math.cos(theta)) / (math.sin(inc) * math.sin(theta)), -1.0, 1.0)
    phi = np.arccos(cphi)
    wgt = (math.pi / n) * 0.5 * (x1 - x2) * np.sin(t) / (math.sin(inc) * math.sin(theta) * np.sqrt(np.maximum(1.0 - cphi * cphi, 1e-300)))
    return np.array([np.sum(wgt * np.cos(m * phi) * x) / math.pi for m in range(mmax + 1)])

# usage: python -I -B -X utf8 count_at_emission_check.py     (a second; no reference material is needed)
# The last line, RESULT ..., is read by 2-python-solution/check.py --slow: the largest difference at 64 nodes over the three cases.
u = 4.135 / 12.0        # L / R for 1.4 Suns and 12 km, the test problems' star
cases = [("never hidden: observer 60 deg, element 20 deg", math.radians(60), math.radians(20)),
         ("hidden for part of the turn: observer 80 deg, element 60 deg", math.radians(80), math.radians(60)),
         ("observer passes straight over it: both 90 deg", math.radians(90), math.radians(90))]
worst = {8: 0.0, 16: 0.0, 32: 0.0, 64: 0.0}
for name, inc, th in cases:
    ref = usual_way(u, inc, th, 6)
    print(name, flush=True)
    print("   usual way, 400,000 phases: " + " ".join("%+.7f" % v for v in ref), flush=True)
    for n in (8, 16, 32, 64):
        got = at_emission(u, inc, th, 6, n)
        worst[n] = max(worst[n], float(np.max(np.abs(got - ref))))
        print("   at emission, %2d nodes:     " % n + " ".join("%+.7f" % v for v in got) + "   largest difference %.1e" % np.max(np.abs(got - ref)), flush=True)
print("RESULT " + " ".join("largest_%d=%.3e" % (n, worst[n]) for n in (8, 16, 32, 64)))
