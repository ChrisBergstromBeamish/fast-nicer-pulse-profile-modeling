# -*- coding: utf-8 -*-
"""Light's path outside a non-spinning mass, from clocks and rulers (22 Nov 2026).

PREMISES (no formula here is taken from general relativity; see NICER_SPEED_PLAN_v0.2.md, section 3):
  P1  clock rate of a frame at rest at radius r, against far clocks: f(r) = sqrt(1 - L/r)     [second paper, Section 3]
  P2  rulers: radial dl = dr / f(r); sideways r dpsi (r is the circumference over 2 pi)        [second paper, Section 3]
  P3  light moves at c in every frame at rest                                                  [second paper, H2]
  P4  along the path b = r sin(alpha) / f(r) stays constant, alpha being the angle between the ray and the
      outward radius as the frame at rest there reads it      [the angular momentum law: not derived in the papers; a named premise, accepted by the author 23 Nov 2026]
CONSEQUENCES (two lines each):
  in the frame at rest, in its time dt_loc: radial step dr/f = c cos(alpha) dt_loc, sideways step r dpsi = c sin(alpha) dt_loc,
  and dt_loc = f dt, so   dpsi/dr = b / (r^2 sqrt(Q)),   c dt/dr = 1 / (f^2 sqrt(Q)),   Q = 1 - b^2 f^2 / r^2 = cos^2(alpha).
UNITS: lengths in units of the star's own radius R at the emission point; u = L/R; x = cos(alpha) at the surface.
  With y = R/r and kb^2 = (b/R)^2 = (1 - x^2)/(1 - u):   Q = P(y) = 1 - kb^2 y^2 (1 - u y).
  swept angle   Theta(u, x) = kb * integral_0^1 dy / sqrt(P)
  delay behind a radial ray from the same radius, in units of R/c:
                tau(u, x)   = kb^2 * integral_0^1 dy / (sqrt(P) (1 + sqrt(P)))
  With w = 1 - y:  P = x^2 + kb^2 w ((2 - 3u) - (1 - 3u) w - u w^2), exactly; this form has no cancellation near the edge.
"""
import math
import numpy as np

_GL = {}


def _gl(n):
    if n not in _GL:
        g, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (g + 1.0), 0.5 * w)
    return _GL[n]


def theta_tau(u, x, n=160):
    """Swept angle Theta and delay tau (units R/c) for outward rays. u = L/R (scalar), x = cos(alpha) in (0, 1] (array).
    Radial variable v = sqrt(1 - R/r); near the star's edge (small x) the variable v = (x/sqrt(a)) sinh(s) keeps the
    integrand smooth, with a = kb^2 (2 - 3u) the slope of P at the surface."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    kb2 = (1.0 - x * x) / (1.0 - u)
    kb = np.sqrt(kb2)
    th = np.empty_like(x)
    ta = np.empty_like(x)
    big = x >= 0.25
    if np.any(big):
        v, w = _gl(n)
        w2 = (v * v)[None, :]
        P = (x[big] ** 2)[:, None] + kb2[big, None] * w2 * ((2.0 - 3.0 * u) - (1.0 - 3.0 * u) * w2 - u * w2 * w2)
        sP = np.sqrt(P)
        th[big] = kb[big] * (((2.0 * v)[None, :] / sP) @ w)
        ta[big] = kb2[big] * (((2.0 * v)[None, :] / (sP * (1.0 + sP))) @ w)
    small = ~big
    if np.any(small):
        g, w = _gl(n)
        xs = np.maximum(x[small], 1e-14)
        a = kb2[small] * (2.0 - 3.0 * u)
        smax = np.arcsinh(np.sqrt(a) / xs)                      # v runs from 0 to 1
        s = smax[:, None] * g[None, :]
        v = (xs / np.sqrt(a))[:, None] * np.sinh(s)
        dv = (xs / np.sqrt(a))[:, None] * np.cosh(s) * smax[:, None]
        w2 = v * v
        P = (xs ** 2)[:, None] + kb2[small, None] * w2 * ((2.0 - 3.0 * u) - (1.0 - 3.0 * u) * w2 - u * w2 * w2)
        sP = np.sqrt(P)
        th[small] = kb[small] * ((2.0 * v * dv / sP) @ w)
        ta[small] = kb2[small] * ((2.0 * v * dv / (sP * (1.0 + sP))) @ w)
    return th, ta


def lensing_factor(u, x, h=1e-5):
    """d cos(alpha) / d cos(Theta), by a centred difference of the forward function (used only by the plain version)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    xp = np.minimum(x + h, 1.0 - 1e-12)
    xm = np.maximum(x - h, 1e-9)
    cp = np.cos(theta_tau(u, xp)[0])
    cm = np.cos(theta_tau(u, xm)[0])
    return (xp - xm) / (cp - cm)

def theta_tau_any(u, x, n=160):
    """Swept angle and delay (units R/c, behind a radial ray from the same radius) for cos(alpha) = x in (-1, 1].
    For x < 0 the ray starts inward, turns at the first radius where Q = 0, and comes back out through R: to the outward
    values at the same b it adds twice the angle and twice the far-clock time of the stretch between the turning radius and R.
    With s = R/r - 1 on that stretch:  P = kb^2 (s_t - s) g(s),  g = (2 - 3u) + (1 - 3u)(s_t + s) - u (s_t^2 + s_t s + s^2),
    which has no cancellation; s = s_t (1 - v^2) removes the square-root end."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    th, ta = theta_tau(u, np.abs(x), n=n)
    neg = x < 0.0
    if np.any(neg):
        xa = -x[neg]
        kb2 = (1.0 - xa * xa) / (1.0 - u)
        c1, c2, c3 = (2.0 - 3.0 * u), (1.0 - 3.0 * u), -u
        st = xa * xa / (kb2 * c1)
        for _ in range(60):                                   # Newton for the turning point: kb^2 s (c1 + c2 s + c3 s^2) = x^2
            F = kb2 * st * (c1 + c2 * st + c3 * st * st) - xa * xa
            dF = kb2 * (c1 + 2.0 * c2 * st + 3.0 * c3 * st * st)
            st = st - F / dF
        v, w = _gl(n)
        s = st[:, None] * (1.0 - v * v)[None, :]
        g = kb2[:, None] * (c1 + c2 * (st[:, None] + s) + c3 * (st[:, None] ** 2 + st[:, None] * s + s * s))
        base = 2.0 * np.sqrt(st)[:, None] / np.sqrt(g)       # ds / sqrt(P) = base dv
        th[neg] += 2.0 * np.sqrt(kb2) * (base @ w)
        ta[neg] += 2.0 * ((base / ((1.0 + s) ** 2 * (1.0 - u * (1.0 + s)))) @ w)
    return th, ta
