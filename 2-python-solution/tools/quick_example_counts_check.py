"""Our expected counts against the expected-counts file that ships with the public code's quick example
(xpsi_good_expected_hreadable.dat: channel, lower phase edge, expected counts; which run of the public code made it is not
stated in the example; its documentation shows a synthesis at the same star with finer numerical settings).
usage: python -I -B -X utf8 quick_example_counts_check.py        (a few seconds; the public example's data from settings.cfg)
The last line, RESULT ..., is read by 2-python-solution/check.py --slow (the converged route's total and largest difference)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quick_example as qe
import public_likelihood as pl
import settings
data_dir = settings.data_dir()
counts, bg, support = qe.load_data(data_dir)
edges_E = pl.example_instrument_edges()
ex = np.loadtxt(os.path.join(data_dir, "xpsi_good_expected_hreadable.dat"))
exp_file = ex[:, 2].reshape(291, 32)
print("file: channels %g to %g, first phase edges %s, total %.3f" % (ex[0, 0], ex[-1, 0], ex[:3, 1], exp_file.sum()))
result = {}
for label, star in (("ours, converged (48 rings per panel, 32 points, 48 energies)", qe.ours(qe.P, counts, support, edges_E, 48, 32, 48, return_star=True)),
                    ("ours, working settings (10 rings per panel, 16 points, 24 energies)", qe.ours(qe.P, counts, support, edges_E, 10, 16, 24, return_star=True)),
                    ("the public code's two interpolation steps on our curve", qe.emulated(qe.P, counts, support, edges_E, 64, 32, 3073, return_star=True))):
    mine = qe.T_EXP / qe.NB * star + bg[:, None] / qe.NB
    src_file = exp_file - bg[:, None] / qe.NB
    src_mine = qe.T_EXP / qe.NB * star
    rel = mine / exp_file - 1.0
    rs = src_mine / src_file - 1.0
    print(label)
    print("   total expected counts: ours %.3f, file %.3f, ratio - 1 = %+.2e" % (mine.sum(), exp_file.sum(), mine.sum() / exp_file.sum() - 1.0))
    print("   per channel and bin, (ours / file) - 1: largest %.2e, r.m.s. %.2e, mean %+.2e ; source part alone: largest %.2e, r.m.s. %.2e"
          % (np.max(np.abs(rel)), np.sqrt(np.mean(rel**2)), rel.mean(), np.max(np.abs(rs)), np.sqrt(np.mean(rs**2))))
    pc = (mine.sum(axis=1) / exp_file.sum(axis=1) - 1.0)
    pb = (mine.sum(axis=0) / exp_file.sum(axis=0) - 1.0)
    print("   per channel (summed over bins): largest %.2e ; per phase bin (summed over channels): largest %.2e" % (np.max(np.abs(pc)), np.max(np.abs(pb))))
    print("   per channel, first 4 and last 4: %s ... %s" % (np.array2string(pc[:4], precision=2), np.array2string(pc[-4:], precision=2)))
    print("   per phase bin: %s" % np.array2string(pb, precision=1, max_line_width=250))
    if not result:
        result = dict(total_ratio=mine.sum() / exp_file.sum() - 1.0, largest_rel=np.max(np.abs(rel)), per_channel_largest=np.max(np.abs(pc)),
                      interior_channels_largest=np.max(np.abs(pc[1:-2])), per_bin_largest=np.max(np.abs(pb)))
print("RESULT total_ratio=%+.3e total_ratio_abs=%.3e largest_rel=%.3e per_channel_largest=%.3e interior_channels_largest=%.3e per_bin_largest=%.3e"
      % (result["total_ratio"], abs(result["total_ratio"]), result["largest_rel"], result["per_channel_largest"], result["interior_channels_largest"], result["per_bin_largest"]))
