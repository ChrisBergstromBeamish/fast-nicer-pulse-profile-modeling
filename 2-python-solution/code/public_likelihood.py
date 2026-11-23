# -*- coding: utf-8 -*-
"""The public code's likelihood, written out again from its source (X-PSI, commit bb5f461; specification:
XPSI_MODEL_SPEC.md sections 8.3, 8.5 and 9). Nothing of X-PSI is imported or run. This file is the definition of the number
to be reproduced; it is not physics of ours and it holds no physics: interpolation, sums, one integral per channel.

  akima_periodic(x, y)          the interpolant of the GSL library called "akima-periodic" (algorithm read from GSL's akima.c)
  energy_integrate(...)         spec 8.3: the signal at 64 energies -> photons per instrument input interval
  phase_integrate(...)          spec 8.5: the count rate at 64 phases -> mean rate in each data phase bin
  marginal_loglike_channel(...) spec 9: ln of the integral over the background rate of the Poisson probability of one
                                channel's counts, with the public code's own search for the limits (reference form)
  marginal_loglike(...)         the sum over channels
"""
import math
import numpy as np


class AkimaPeriodic(object):
    """GSL's gsl_interp_akima_periodic through (x[k], y[..., k]). The end slopes are wrapped around whether or not the data
    are periodic (the public code uses this type for the energy integration too, where the data are not periodic)."""

    def __init__(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        N = x.size
        h = np.diff(x)
        m = np.diff(y, axis=-1) / h                                   # m[0 .. N-2]
        me = np.concatenate([m[..., N - 3:N - 1], m, m[..., 0:2]], axis=-1)   # me[k + 2] = m[k], k = -2 .. N
        i = np.arange(N - 1) + 2
        d_next = np.abs(me[..., i + 1] - me[..., i])                   # |m[i+1] - m[i]|
        d_prev = np.abs(me[..., i - 1] - me[..., i - 2])               # |m[i-1] - m[i-2]|
        NE = d_next + d_prev
        d_nn = np.abs(me[..., i + 2] - me[..., i + 1])                 # |m[i+2] - m[i+1]|
        d_pp = np.abs(me[..., i] - me[..., i - 1])                     # |m[i] - m[i-1]|
        NE_next = d_nn + d_pp
        mi, mim1, mip1 = me[..., i], me[..., i - 1], me[..., i + 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha_i = np.where(NE != 0.0, d_prev / np.where(NE != 0.0, NE, 1.0), 0.0)
            alpha_ip1 = np.where(NE_next != 0.0, d_pp / np.where(NE_next != 0.0, NE_next, 1.0), 0.0)
        tL = np.where(NE_next != 0.0, (1.0 - alpha_ip1) * mi + alpha_ip1 * mip1, mi)
        b = (1.0 - alpha_i) * mim1 + alpha_i * mi
        c = (3.0 * mi - 2.0 * b - tL) / h
        d = (b + tL - 2.0 * mi) / (h * h)
        flat = NE == 0.0
        self.b = np.where(flat, mi, b)
        self.c = np.where(flat, 0.0, c)
        self.d = np.where(flat, 0.0, d)
        self.x, self.y, self.h = x, y, h
        # primitive at the nodes, from x[0]
        seg = h * (y[..., :-1] + h * (0.5 * self.b + h * (self.c / 3.0 + 0.25 * self.d * h)))
        self.P = np.concatenate([np.zeros(y.shape[:-1] + (1,)), np.cumsum(seg, axis=-1)], axis=-1)

    def _index(self, q):
        return np.clip(np.searchsorted(self.x, q, side="right") - 1, 0, self.x.size - 2)

    def __call__(self, q):
        q = np.asarray(q, dtype=float)
        k = self._index(q)
        t = q - self.x[k]
        return self.y[..., k] + t * (self.b[..., k] + t * (self.c[..., k] + t * self.d[..., k]))

    def primitive(self, q):
        """Integral of the interpolant from x[0] to q (q: array; the result has shape y.shape[:-1] + q.shape)."""
        q = np.asarray(q, dtype=float)
        k = self._index(q)
        t = q - self.x[k]
        return self.P[..., k] + t * (self.y[..., k] + t * (0.5 * self.b[..., k] + t * (self.c[..., k] / 3.0
                                                                                      + 0.25 * self.d[..., k] * t)))

    def integral(self, a, b):
        return self.primitive(b) - self.primitive(a)


def example_instrument_edges():
    """The 292 edges of the quick example's 291 instrument input intervals, made the way its main.py makes them
    (repeated addition of 0.005 from 0.1 and 0.105; rows 10 to 300 kept)."""
    low, high = [0.1], [0.105]
    while low[-1] <= 14.995:
        low.append(low[-1] + 0.005)
        high.append(high[-1] + 0.005)
    return np.array([low[10]] + high[10:301])


def example_energies(edges, num=64):
    return np.logspace(np.log10(edges[0]), np.log10(edges[-1]), num, base=10.0)


def energy_integrate(signal, energies, edges):
    """signal[e, k]: photons/cm^2/s/keV at energies[e] and phase k. Returns binned[q, k]: photons/cm^2/s in the input interval
    [edges[q], edges[q+1]], by the exact integral of the akima-periodic interpolant of E F ln(10) in log10(E)."""
    x = np.log10(energies)
    le = np.log10(edges)
    y = (np.power(10.0, x)[:, None] * signal * math.log(10.0)).T          # (phase, energy)
    ak = AkimaPeriodic(x, y)
    up = np.minimum(le[1:], x[-1])
    return (ak.primitive(up) - ak.primitive(le[:-1])).T                    # (interval, phase)


def phase_integrate(rate, signal_phases, bin_edges, phase_shift=0.0):
    """rate[r, k]: counts/s in channel r at signal_phases[k] (cycles, 0 to 1 inclusive). Returns star[r, j] = n times the
    integral over data bin j of the akima-periodic interpolant (the mean rate in the bin), with the public code's rules:
    a piece is added only if positive, and a negative total is set to zero."""
    ak = AkimaPeriodic(signal_phases, rate)
    nb = len(bin_edges) - 1
    out = np.zeros(rate.shape[:-1] + (nb,))
    for j in range(nb):
        pa = bin_edges[j] + phase_shift
        pb = bin_edges[j + 1] + phase_shift
        if abs(pb - pa - 1.0) < 1e-12:
            pa, pb = 0.0, 1.0
        else:
            pa -= math.floor(pa)
            pb -= math.floor(pb)
        if pa < pb:
            v = ak.integral(pa, pb)
            out[..., j] += np.where(v > 0.0, v, 0.0)
        else:
            v = ak.integral(pa, 1.0)
            out[..., j] += np.where(v > 0.0, v, 0.0)
            v = ak.integral(0.0, pb)
            out[..., j] += np.where(v > 0.0, v, 0.0)
    return np.maximum(out, 0.0) * nb


_GL = {}


def _gl(n):
    if n not in _GL:
        g, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (g + 1.0), 0.5 * w)
    return _GL[n]


def _log_cond(B, d, s, SCALE):
    """sum_j [d_j ln(c_j) - c_j], c_j = SCALE (s_j + B), for an array of B; -inf where the public code's integrand is zero."""
    c = SCALE * (s[None, :] + np.asarray(B, dtype=float)[:, None])
    pos = c > 0.0
    both0 = (np.abs(c) < 1e-12) & (np.abs(d[None, :]) < 1e-12)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(pos, d[None, :] * np.log(np.where(pos, c, 1.0)) - c, 0.0)
    bad = ~(pos | both0)
    out = term.sum(axis=1)
    out[bad.any(axis=1)] = -np.inf
    return out


def channel_limits(d, s, T, support, epsilon=1.0e-3, sigmas=10.0):
    """The public code's search for the best background rate and the limits of the integral (spec 9.5, 9.6).
    Returns (lower, upper, B_for_integrand, B_best)."""
    n = d.size
    SCALE = T / n
    av_STAR = s.sum() / n
    av_DATA = d.sum() / T
    has = d > 0.0

    def dlog(B):
        if np.any(has & (s + B <= 0.0)):
            return math.inf
        return float(np.sum(d[has] / (s[has] + B))) - T

    def delta(B):
        y = float(np.sum(d[has] / (s[has] + B)))
        x = float(np.sum(d[has] / (s[has] + B) ** 2))
        y = 2.0 * T - 2.0 * y
        x *= 2.0
        return -y / x, math.sqrt(2.0 / x)

    B_min = 0.0
    B = av_DATA - av_STAR
    if B <= B_min:
        min_counts = -1.0
        if np.any(np.abs(s) < 1e-12):
            min_counts = -2.0
        if min_counts == -2.0:
            for j in range(n):
                if (min_counts == -2.0 and d[j] > 0.0) or (0.0 < d[j] < min_counts):
                    min_counts = d[j]
        if min_counts != -1.0:
            if min_counts == -2.0:
                B, B_min = 0.0, 0.0
            else:
                B = 0.01 * min_counts / SCALE
                B_min = 0.1 * B
        else:
            B = B_min
    g = dlog(B_min)
    if g <= 0.0:
        B = B_min
    else:
        B_low, B_high = B_min, av_DATA
        for j in range(n):
            if s[j] < 0.0 and av_DATA - s[j] > B_high:
                B_high = av_DATA - s[j]
        if B < B_low or B >= B_high:
            B = 0.5 * (B_low + B_high)
        dB, std = delta(B)
        counter = 0
        while not abs(dB) <= epsilon * std and counter < 1000:
            g = dlog(B)
            if g > 0.0:
                B_low = B
            elif g < 0.0:
                B_high = B
            else:
                break
            B_new = B + dB
            if not (B_low < B_new < B_high):
                B_new = 0.5 * (B_low + B_high)
            B = B_new
            dB, std = delta(B)
            counter += 1
    se = float(np.sum(d[has] / (s[has] + B) ** 2))
    std_est = math.sqrt(1.0 / se) if se > 0.0 else 1e90
    lower = B - sigmas * std_est
    upper = B + sigmas * std_est
    if lower < B_min:
        lower = B_min
    B_for = B
    if lower < support[0]:
        lower = support[0]
        if upper < support[0] and support[1] > 0.0:
            upper = support[1]
            B_for = support[0]
        elif upper < support[0]:
            upper = support[0] + sigmas * std_est
            B_for = support[0]
    if upper > support[1] and support[1] > 0.0:
        upper = support[1]
        if lower > support[1]:
            lower = support[0]
            B_for = support[1]
    if B_for < lower:
        B_for = lower
    elif B_for > upper:
        B_for = upper
    if dlog(upper) >= 0.0:
        B_for = upper
    elif dlog(lower) <= 0.0:
        B_for = lower
    return lower, upper, B_for, B


def marginal_loglike_channel(d, s, T, support, epsilon=1.0e-3, sigmas=10.0, rtol=1.0e-12):
    """ln of INT_lower^upper PROD_j Poisson(d_j | (T/n)(s_j + B)) dB for one channel (reference form: the public code's
    limits; the integral by Gauss-Legendre panels, doubled until two answers agree to rtol)."""
    d = np.asarray(d, dtype=float)
    s = np.asarray(s, dtype=float)
    n = d.size
    SCALE = T / n
    fact = -float(sum(math.lgamma(v + 1.0) for v in d))
    if abs(d.sum() / T) < 1e-12 and abs(s.sum() / n) < 1e-12:
        lower = max(0.0, support[0])
        upper = 10.0 / T
        if upper > support[1] and support[1] > 0.0:
            upper = support[1]
        return math.log((math.exp(-lower * T) - math.exp(-upper * T)) / T) + fact
    lower, upper, B_for, _ = channel_limits(d, s, T, support, epsilon, sigmas)
    if d.sum() == 0.0:
        # no counts in the channel: the integrand is exp(-(T/n) sum s - T B) exactly
        return -SCALE * float(s.sum()) - T * lower + math.log(-math.expm1(-T * (upper - lower)) / T) + fact
    A = float(_log_cond(np.array([B_for]), d, s, SCALE)[0])
    if not math.isfinite(A):
        return -math.inf                                   # the public code's failure exit (a number near -1e90)
    # The logarithm of the integrand is concave in B and largest at B_for. Only the part of [lower, upper] within 60 of that
    # largest value is integrated; what is left out is below e^-60 of the largest value times its own width.
    has = d > 0.0
    c2 = float(np.sum(d[has] / (s[has] + B_for) ** 2))
    g1 = abs(float(np.sum(d[has] / (s[has] + B_for))) - T)
    step0 = min(v for v in (1.0 / math.sqrt(c2) if c2 > 0.0 and math.isfinite(c2) else math.inf,
                            1.0 / g1 if g1 > 0.0 and math.isfinite(g1) else math.inf, upper - lower) if v > 0.0)
    lo, hi = B_for, B_for
    step = step0
    while hi < upper:
        hi = min(upper, hi + step)
        if A - float(_log_cond(np.array([hi]), d, s, SCALE)[0]) >= 60.0:
            break
        step *= 2.0
    step = step0
    while lo > lower:
        lo = max(lower, lo - step)
        if A - float(_log_cond(np.array([lo]), d, s, SCALE)[0]) >= 60.0:
            break
        step *= 2.0
    g, w = _gl(20)
    prev = None
    panels = 16
    while True:
        ed = np.linspace(lo, hi, panels + 1)
        Bn = (ed[:-1, None] + (ed[1:] - ed[:-1])[:, None] * g[None, :]).ravel()
        wn = ((ed[1:] - ed[:-1])[:, None] * w[None, :]).ravel()
        val = float(np.sum(wn * np.exp(_log_cond(Bn, d, s, SCALE) - A)))
        if prev is not None and abs(val - prev) <= rtol * abs(val):
            break
        if panels >= 16384:
            raise ArithmeticError("the reference integral did not settle")
        prev = val
        panels *= 2
    if not (val > 0.0 and math.isfinite(val)):
        raise ArithmeticError("the reference integral came out zero or not finite")
    return math.log(val) + A + fact


def marginal_loglike(counts, star, T, support, epsilon=1.0e-3, sigmas=10.0):
    """The public code's returned log-likelihood: the sum over channels (spec 9.7). Also the slim exit's test value."""
    counts = np.asarray(counts, dtype=float)
    n = counts.shape[1]
    tot = 0.0
    for r in range(counts.shape[0]):
        src = star[r].sum() * T / n
        if src - 20.0 * math.sqrt(src) - counts[r].sum() > 0.0:
            raise ValueError("channel %d: the public code would leave by its early exit (model far above the data)" % r)
        tot += marginal_loglike_channel(counts[r], star[r], T, support[r], epsilon, sigmas)
    return tot


# ---------------------------------------------------------------------------------------------------------------------
# The same number, computed for all channels at once (the fast form).
# For each channel: the best background rate by Newton's steps (all channels together), the public code's limits, then one
# Gauss-Legendre rule on the part of [lower, upper] where the integrand is above e^-36 of its largest value. The number of
# points is chosen from how far the logarithm of the integrand falls over that part: 10 points to a fall of 0.5, 14 to 2,
# 20 to 8, 28 to 18, 40 to 45, 64 beyond (measured against the step-by-step form on 247 hard made-up channels: 6e-11 per
# channel; a measurement on those channels, not a bound for every channel).
# Channels with a bin of zero model rate, or with no counts, go to the reference form above.
def log_factorial_sums(counts):
    """-sum_j ln(d_j!) per channel. Depends on the data only: computed once per data set."""
    c = np.asarray(counts)
    table = np.array([math.lgamma(k + 1.0) for k in range(int(c.max()) + 1)])
    return -table[c.astype(np.int64)].sum(axis=1)


_NODE_CLASSES = ((0.5, 10), (2.0, 14), (8.0, 20), (18.0, 28), (45.0, 40), (1e300, 64))
DROP = 36.0
CUT_FALL = 32.0


def marginal_loglike_fast(counts, star, T, support, fact, sigmas=10.0, details=None, epsilon=1.0e-3):
    d = np.asarray(counts, dtype=float)
    s = np.asarray(star, dtype=float)
    nch, n = d.shape
    SCALE = T / n
    dsum = d.sum(axis=1)
    src = s.sum(axis=1) * SCALE
    if np.any(src - 20.0 * np.sqrt(np.maximum(src, 0.0)) - dsum > 0.0):
        return -math.inf                                   # the public code leaves here with a number near -1e90
    usual = np.all(s > 1e-12, axis=1) & (dsum > 0.0)       # the public code treats a rate below 1e-12 as zero
    out = np.empty(nch)
    for r in np.where(~usual)[0]:
        out[r] = marginal_loglike_channel(d[r], s[r], T, support[r], epsilon=epsilon, sigmas=sigmas)
    iu = np.where(usual)[0]
    if iu.size == 0:
        return float(out.sum())
    d, s, sup = d[iu], s[iu], np.asarray(support, dtype=float)[iu]
    # The best background rate, by the public code's own steps (spec 9.5), for all channels together: its start, its
    # bracket, its Newton step or a bisection, and its stop at |step| <= epsilon * sigma. The iterate it stops at is kept,
    # because the limits B +- sigmas * sigma are built from it (where a limit cuts through a part of the integrand that
    # matters, the returned number depends on that iterate at the 1e-6 level).
    av_data = dsum[iu] / T
    B = av_data - s.sum(axis=1) / n
    B = np.where(B <= 0.0, 0.0, B)
    r0 = d / s
    climbing = (r0.sum(axis=1) - T) > 0.0                           # the slope at B = 0; if not positive the best rate is 0
    B_low = np.zeros(iu.size)
    B_high = av_data.copy()
    B = np.where(climbing & ((B < B_low) | (B >= B_high)), 0.5 * (B_low + B_high), B)
    B = np.where(climbing, B, 0.0)
    r1 = d / (s + B[:, None])
    y = r1.sum(axis=1)
    x2 = (r1 / (s + B[:, None])).sum(axis=1)
    dB = (y - T) / x2
    active = climbing & ~(np.abs(dB) <= epsilon / np.sqrt(x2))
    for _ in range(1000):
        if not np.any(active):
            break
        g_ = y - T
        active = active & (g_ != 0.0)
        B_low = np.where(active & (g_ > 0.0), B, B_low)
        B_high = np.where(active & (g_ < 0.0), B, B_high)
        B_new = B + dB
        B_new = np.where((B_low < B_new) & (B_new < B_high), B_new, 0.5 * (B_low + B_high))
        B = np.where(active, B_new, B)
        r1 = d / (s + B[:, None])
        y = r1.sum(axis=1)
        x2 = (r1 / (s + B[:, None])).sum(axis=1)
        dB = (y - T) / x2
        active = active & ~(np.abs(dB) <= epsilon / np.sqrt(x2))
    std = 1.0 / np.sqrt(x2)
    # the public code's limits (spec 9.5), in its order
    lower = np.maximum(B - sigmas * std, 0.0)
    upper = B + sigmas * std
    s0, s1 = sup[:, 0], sup[:, 1]
    has_up = s1 > 0.0
    below = lower < s0
    out_low = below & (upper < s0)
    upper = np.where(out_low & has_up, s1, np.where(out_low, s0 + sigmas * std, upper))
    lower = np.where(below, s0, lower)
    over = has_up & (upper > s1)
    out_high = over & (lower > s1)
    upper = np.where(over, s1, upper)
    lower = np.where(out_high, s0, lower)
    # the largest value of the integrand on [lower, upper], and the part that matters
    Bs = np.clip(B, lower, upper)
    rs = d / (s + Bs[:, None])
    gs = rs.sum(axis=1) - T
    cs = (rs / (s + Bs[:, None])).sum(axis=1)
    root = np.sqrt(gs * gs + 2.0 * DROP * cs)
    # where g x - c x^2/2 = -DROP, to the right and to the left; each written in the form without cancellation
    xR = np.where(gs >= 0.0, (root + gs) / cs, 2.0 * DROP / (root - gs))
    xL = np.where(gs <= 0.0, (root - gs) / cs, 2.0 * DROP / (root + gs))
    lo = np.maximum(lower, Bs - xL)
    hi = np.minimum(upper, Bs + xR)
    ssum = SCALE * s.sum(axis=1)
    ls = (d * np.log(SCALE * (s + Bs[:, None]))).sum(axis=1) - T * Bs - ssum

    def lfun(idx, Bq):
        return (d[idx] * np.log(SCALE * (s[idx] + Bq[:, None]))).sum(axis=1) - T * Bq - ssum[idx]
    # The two cut points came from a parabola through the largest value. The true logarithm is checked there, and a cut is
    # moved outward until the true fall is at least CUT_FALL (the logarithm is concave, so beyond the cut it only falls
    # faster: what is left out is below e^-CUT_FALL of the largest value times a width of the order of the kept part).
    fall_hi = np.zeros(iu.size)
    fall_lo = np.zeros(iu.size)
    for _ in range(60):
        idx = np.where((hi < upper) & (fall_hi < CUT_FALL))[0]
        if idx.size == 0:
            break
        fall_hi[idx] = ls[idx] - lfun(idx, hi[idx])
        short = idx[fall_hi[idx] < CUT_FALL]
        hi[short] = np.minimum(upper[short], Bs[short] + 1.5 * (hi[short] - Bs[short]))
    for _ in range(60):
        idx = np.where((lo > lower) & (fall_lo < CUT_FALL))[0]
        if idx.size == 0:
            break
        fall_lo[idx] = ls[idx] - lfun(idx, lo[idx])
        short = idx[fall_lo[idx] < CUT_FALL]
        lo[short] = np.maximum(lower[short], Bs[short] - 1.5 * (Bs[short] - lo[short]))
    span = np.maximum(hi - Bs, Bs - lo)
    fall = np.maximum(np.abs(gs) * span + 0.5 * cs * span * span, np.maximum(fall_hi, fall_lo))
    res = np.empty(iu.size)
    # Where the part that matters is narrow against every bin's rate, the logarithms are summed as a power series about Bs:
    #   l(Bs + x) - l(Bs) = (M_1 - T) x - M_2 x^2/2 + M_3 x^3/3 - ...,   M_k = sum_j d_j / (s_j + Bs)^k,
    # cut after 8 terms; the cut-off is below  D rho^9 / (9 (1 - rho)),  D = the channel's counts, rho = largest |x| / (s_j + Bs).
    rinv = 1.0 / (s + Bs[:, None])
    rho = span * rinv.max(axis=1)
    series = (dsum[iu] * rho**9 / (9.0 * (1.0 - np.minimum(rho, 0.5))) < 1e-10) & (fall <= 2.0)
    if np.any(series):
        ss = np.where(series)[0]
        rr = rinv[ss]
        pw = d[ss] * rr
        coef = np.empty((8, ss.size))
        coef[0] = pw.sum(axis=1) - T
        for k_ in range(2, 9):
            pw = pw * rr
            coef[k_ - 1] = pw.sum(axis=1) * ((-1.0) ** (k_ + 1) / k_)
        g, w = _gl(10)
        x = (lo[ss] - Bs[ss])[:, None] + (hi[ss] - lo[ss])[:, None] * g[None, :]
        acc = coef[7][:, None] * x
        for k_ in range(6, -1, -1):
            acc = (acc + coef[k_][:, None]) * x
        res[ss] = np.log(((hi[ss] - lo[ss])[:, None] * w[None, :] * np.exp(acc)).sum(axis=1)) + ls[ss]
    fall = np.where(series, -1.0, fall)
    prev = 0.0
    for limit, npts in _NODE_CLASSES:
        sel = np.where((fall > prev) & (fall <= limit))[0] if prev > 0.0 else np.where((fall >= 0.0) & (fall <= limit))[0]
        prev = limit
        if sel.size == 0:
            continue
        g, w = _gl(npts)
        Bn = lo[sel, None] + (hi[sel] - lo[sel])[:, None] * g[None, :]                 # (m, npts)
        lg = np.log(SCALE * (s[sel][:, None, :] + Bn[:, :, None]))                       # (m, npts, n)
        ll = (lg * d[sel][:, None, :]).sum(axis=2) - T * Bn - (SCALE * s[sel].sum(axis=1))[:, None]
        res[sel] = np.log(((hi[sel] - lo[sel])[:, None] * w[None, :] * np.exp(ll - ls[sel, None])).sum(axis=1)) + ls[sel]
    if not np.all(np.isfinite(res)):
        raise ArithmeticError("a channel's integral came out zero or not finite in the fast form")
    out[iu] = res + np.asarray(fact)[iu]
    if details is not None:
        details.update(B=B, lower=lower, upper=upper, lo=lo, hi=hi, fall=fall, usual=usual, per_channel=out, series=series)
    return float(out.sum())
