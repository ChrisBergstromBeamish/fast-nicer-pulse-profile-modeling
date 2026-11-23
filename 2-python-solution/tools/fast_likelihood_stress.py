"""The all-channels-at-once likelihood against the reference form, channel by channel, on made-up channels chosen to be hard:
few counts, source far above or below the data, wide and narrow limits on the background, Codex's round-3 counterexample
(32 equal source rates of 0.0005 per second, one count, no upper limit).
usage: python -I -B -X utf8 fast_likelihood_stress.py        (a few seconds; no reference material is needed)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import public_likelihood as pl
rng = np.random.default_rng(7)
T, n = 1000.0, 32
rows_d, rows_s, rows_sup, labels = [], [], [], []
# Codex's counterexample
d = np.zeros(n); d[5] = 1.0
rows_d.append(d); rows_s.append(np.full(n, 0.0005)); rows_sup.append([0.0, -1.0]); labels.append("counterexample")
for level in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0):          # mean counts per bin
    for ratio in (0.05, 0.3, 0.8, 1.0, 1.2, 3.0):                             # source over data
        for sup in ("none", "narrow", "wide", "low", "high"):
            prof = 1.0 + 0.9 * np.cos(2 * np.pi * (np.arange(n) + 0.5) / n + rng.uniform(0, 6.28))
            lam_d = level * prof
            d = rng.poisson(lam_d).astype(float)
            s = ratio * lam_d / (T / n) * rng.uniform(0.9, 1.1)
            bmode = max(d.sum() / T - s.mean(), 0.0)
            scale = max(bmode, math.sqrt(max(d.sum(), 1.0)) / T)
            S = {"none": [0.0, -1.0], "narrow": [0.9 * scale, 1.1 * scale], "wide": [0.1 * scale, 10.0 * scale],
                 "low": [0.0, 0.2 * scale], "high": [3.0 * scale, 6.0 * scale]}[sup]
            rows_d.append(d); rows_s.append(s); rows_sup.append(S); labels.append("%g per bin, source/data %.2f, limits %s" % (level, ratio, sup))
D, Sr, Sup = np.array(rows_d), np.array(rows_s), np.array(rows_sup)
src = Sr.sum(axis=1) * T / n
screen = src - 20.0 * np.sqrt(src) - D.sum(axis=1) > 0.0          # the public code's early exit: such a channel ends the evaluation
print("%d made-up channels; %d of them would trigger the public code's early exit and are left out" % (D.shape[0], int(screen.sum())))
D, Sr, Sup = D[~screen], Sr[~screen], Sup[~screen]
labels = [l for l, k in zip(labels, screen) if not k]
fact = pl.log_factorial_sums(D)
det = {}
tot = pl.marginal_loglike_fast(D, Sr, T, Sup, fact, details=det)
if not np.isfinite(tot):
    print("the fast form left by the early exit; removing those channels")
ref = np.array([pl.marginal_loglike_channel(D[r], Sr[r], T, Sup[r]) for r in range(D.shape[0])])
diff = det["per_channel"] - ref
bad = np.where(~(np.abs(diff) < 1e-8))[0]
largest = np.nanmax(np.abs(diff))
print("%d made-up channels: largest difference fast minus reference %.2e ; channels above 1e-8: %d ; sent to the reference form by the fast form itself: %d"
      % (D.shape[0], largest, bad.size, int((~det["usual"]).sum())))
for r in bad[:25]:
    print("   %-55s counts %5.0f: fast %.10f, reference %.10f, difference %+.2e" % (labels[r], D[r].sum(), det["per_channel"][r], ref[r], diff[r]))
print("RESULT channels=%d largest=%.3e above=%d" % (D.shape[0], largest if np.all(np.isfinite(diff)) else math.inf, bad.size))
