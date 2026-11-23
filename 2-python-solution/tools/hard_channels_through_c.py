"""The 247 hard made-up channels of fast_likelihood_stress.py (few counts, source far above or below the data, wide, narrow
and displaced limits, Codex's one-count counterexample) through the compiled copy's likelihood, channel by channel, against
the step-by-step Python form. The channels are written to 3-compiled-solution/hard_channels.bin, the compiled program writes
its answers next to it, and both files are removed.
usage: python -I -B -X utf8 hard_channels_through_c.py [nicer_fast | nicer_fast_strict]     (nicer_fast if no name is given)
The last line, RESULT ..., is read by 3-compiled-solution/check.py."""
import math, os, struct, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import public_likelihood as pl, settings
program = sys.argv[1] if len(sys.argv) > 1 else "nicer_fast"
rng = np.random.default_rng(7)
T, n = 1000.0, 32
rows_d, rows_s, rows_sup, labels = [], [], [], []
d = np.zeros(n); d[5] = 1.0
rows_d.append(d); rows_s.append(np.full(n, 0.0005)); rows_sup.append([0.0, -1.0]); labels.append("counterexample")
for level in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0):
    for ratio in (0.05, 0.3, 0.8, 1.0, 1.2, 3.0):
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
keep = ~(src - 20.0 * np.sqrt(src) - D.sum(axis=1) > 0.0)
D, Sr, Sup = D[keep], Sr[keep], Sup[keep]
labels = [l for l, k in zip(labels, keep) if k]
fact = pl.log_factorial_sums(D)
cdir = settings.compiled_dir()
path = os.path.join(cdir, "hard_channels.bin")
with open(path, "wb") as fh:
    fh.write(struct.pack("<iid", D.shape[0], n, T))
    for a in (D, Sr, Sup, fact):
        fh.write(np.ascontiguousarray(a, dtype=np.float64).tobytes())
out = settings.run_compiled(program, ["tables.bin", "channels", "hard_channels.bin"])
print("the compiled copy said:", out.strip() or "(nothing: is it built?)")
if not os.path.isfile(path + ".out"):
    os.remove(path)
    print("RESULT channels=%d largest=inf above=%d largest_vs_fast=inf" % (D.shape[0], D.shape[0]))
    sys.exit(1)
cvals = np.fromfile(path + ".out", dtype=np.float64)
ref = np.array([pl.marginal_loglike_channel(D[r], Sr[r], T, Sup[r]) for r in range(D.shape[0])])
det = {}
pl.marginal_loglike_fast(D, Sr, T, Sup, fact, details=det)
diff = cvals - ref
bad = np.where(~(np.abs(diff) < 1e-8))[0]
largest, largest_fast = np.nanmax(np.abs(diff)), np.nanmax(np.abs(cvals - det["per_channel"]))
print("%d hard channels: compiled minus the step-by-step Python form: largest %.2e ; channels above 1e-8: %d ; compiled minus the fast Python form: largest %.2e"
      % (D.shape[0], largest, bad.size, largest_fast))
for r in bad[:12]:
    print("   %-52s counts %5.0f: compiled %.10f, reference %.10f, difference %+.2e" % (labels[r], D[r].sum(), cvals[r], ref[r], diff[r]))
os.remove(path); os.remove(path + ".out")
print("RESULT channels=%d largest=%.3e above=%d largest_vs_fast=%.3e" % (D.shape[0], largest if cvals.size == D.shape[0] else math.inf, bad.size, largest_fast))
