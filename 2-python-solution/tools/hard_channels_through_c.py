"""The 247 hard made-up channels of fast_likelihood_stress.py (few counts, source far above or below the data, wide, narrow
and displaced limits, Codex's one-count counterexample) through the compiled copy's likelihood, channel by channel, against
the step-by-step Python form."""
import math, os, struct, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import public_likelihood as pl
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
cdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")
path = os.path.join(cdir, "hard_channels.bin")
with open(path, "wb") as fh:
    fh.write(struct.pack("<iid", D.shape[0], n, T))
    for a in (D, Sr, Sup, fact):
        fh.write(np.ascontiguousarray(a, dtype=np.float64).tobytes())
wsl_c = (lambda p: "/mnt/" + p[0].lower() + p[2:].replace(chr(92), "/") if p[1:2] == ":" else p)(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "3-compiled-solution")))  # the same folder as the Linux side of a Windows machine sees it
out = subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", "cd %s && ./nicer_fast tables.bin channels hard_channels.bin" % wsl_c],
                     capture_output=True, text=True).stdout.replace(chr(0), "")
print("the compiled copy said:", out.strip())
cvals = np.fromfile(path + ".out", dtype=np.float64)
ref = np.array([pl.marginal_loglike_channel(D[r], Sr[r], T, Sup[r]) for r in range(D.shape[0])])
det = {}
pl.marginal_loglike_fast(D, Sr, T, Sup, fact, details=det)
diff = cvals - ref
bad = np.where(~(np.abs(diff) < 1e-8))[0]
print("%d hard channels: compiled minus the step-by-step Python form: largest %.2e ; channels above 1e-8: %d ; compiled minus the fast Python form: largest %.2e"
      % (D.shape[0], np.nanmax(np.abs(diff)), bad.size, np.nanmax(np.abs(cvals - det["per_channel"]))))
for r in bad[:12]:
    print("   %-52s counts %5.0f: compiled %.10f, reference %.10f, difference %+.2e" % (labels[r], D[r].sum(), cvals[r], ref[r], diff[r]))
os.remove(path); os.remove(path + ".out")
