# What the comparison with the public code printed on 23 November 2026

The four files in this folder are what the comparison tools printed on 23 November 2026, kept byte for byte as printed
(three of them end their lines with CR LF, as they were written on Windows; `.gitattributes` keeps them so). Every number
the pages take from them is a row of `../CLAIMS.md`, group H.

The public code X-PSI (version 3.5.0b1, at the commit the setup page names) was installed in its own folder on the Linux
side of the same laptop (Intel Core Ultra 7 258V) and run there, one processor core, one thread. Ours ran on the same
laptop, one core, one thread: in Python on the Windows side, compiled on the Linux side. The files were printed by the
project's working copies of the tools named below, which differ from the committed ones only in how they find folders
and take their arguments.

`PY` stands for `python -I -B -X utf8`, run from the repository's top. "Needs X-PSI" means that X-PSI must be installed and
run first (`../README.md`, section 6).

| file | what it is | power state | the command that prints it |
|---|---|---|---|
| `compare_with_public_code_2026-11-23_mains.txt` | X-PSI at its shipped settings and at three finer ones, for the quick example and for its second shipped check, with ours beside it: values, bins, X-PSI's count rates through our likelihood, and for the quick example the times and their ratios. The times to lead with. | mains power (its first line: Windows' battery status 2, "on mains"); X-PSI's times were taken three minutes before ours, in the same state | needs X-PSI: `bash 2-python-solution/tools/run_public_code_timing.sh` on the Linux side, then `PY 2-python-solution/tools/compare_with_public_code.py 2-python-solution/results/public_code_times.jsonl "<the power state>"` |
| `compare_with_public_code_2026-11-23_battery.txt` | The same comparison, taken earlier the same day. Its values and bins are the same as on mains power; its times are longer. | battery power (its first line: battery status 1); X-PSI's times were taken about ten minutes before ours, in the same state | the same as the file above |
| `compare_refine_2026-11-23.txt` | X-PSI with its settings refined one at a time from its shipped ones (its energies, its rotation steps, its mesh, its rays) and at two finer combinations, for both checks, with ours beside it: values and bins. The times in it are from one evaluation each and are not used. | battery power (its first line says so) | needs X-PSI: `bash 2-python-solution/tools/refine_public_code.sh` on the Linux side, then `PY 2-python-solution/tools/compare_with_public_code.py 2-python-solution/results/public_code_refine.jsonl "on battery; values and bins only are used from this file"` |
| `ours_convergence_both_stars_2026-11-23.txt` | Ours at settings finer than ours converged (48 rings per panel, 32 points, 48 energies), at the quick example's star and at the second check's: ln L and the bins against the converged setting. Only our programs run. | not recorded; it holds no time | `PY 2-python-solution/tools/ours_convergence_both_stars.py` (does not need X-PSI) |
