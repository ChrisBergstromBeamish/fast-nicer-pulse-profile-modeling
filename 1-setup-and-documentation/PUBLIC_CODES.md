# The public codes this work is set beside

Two public codes solve the same kind of problem as the programs in this repository: X-PSI, the public code this work is
compared with, and the graphics-card code of Zhou and Huang. This page says what each is, what of it was used here and what
was not, what each says of its own speed, what was measured of X-PSI on the laptop our programs run on, and why the
comparison with the graphics-card code is not like for like. NASA, the
NICER team and the authors of these codes have not asked for, reviewed or endorsed this work.

## 1. X-PSI

### 1.1 What it is

Its `README.rst`, at the commit read (`bb5f461e605b22c6251b155b27a914a558757f93`), opens: "An open-source package for
neutron star X-ray Pulse Simulation and Inference." And then: "X-PSI is designed to simulate rotationally-modified (pulsed)
surface X-ray emission from neutron stars, taking into account relativistic effects on the emitted radiation. This can then
be used to perform Bayesian statistical inference on real or simulated astronomical data sets." Its source is at the
address its `pyproject.toml` and `README.rst` print, `https://github.com/xpsi-group/xpsi`. Its licence is in
`LICENCE_NOTES.md`.

### 1.2 What of it is used here, and what is not

Used, as targets and as inputs of the comparison:

- the data files of its shipped quick example (`examples/examples_fast/Data`), read as numbers at run time;
- the example's inputs: its star, its instrument, its exposure, its constants, including the flattened shape of its star,
  which `COVERAGE.md` (section 1.4) flags as an input of the example;
- its shipped expected values: the quick example's expected likelihood and tolerance; the expected likelihood and
  tolerance of its second test (`xpsi/tests/test_likelihood.py`); two unit-test values of its likelihood step
  (`xpsi/tests/test_background_marginalisation.py`); and the example's file of expected counts;
- its likelihood definition, written anew by us from its source (`2-python-solution/code/public_likelihood.py`): the
  number to be reproduced, not physics of ours.

Not used by our programs: its code. Nothing of it is imported by our programs or by the checks, and no file of it is
copied here (`LICENCE_NOTES.md`). What the agreement with its shipped values does and does not show is in `COVERAGE.md`,
section 5.

Installed and run beside our programs, for the comparison only: on 23 November 2026, on the author's decision, X-PSI
version 3.5.0b1, at the commit named above, was installed in its own folder on the Linux side of the same laptop by
`2-python-solution/tools/install_reference_code.sh` (X-PSI's documented steps), and run there by
`run_public_code_timing.sh` and `refine_public_code.sh` in the same folder. Both call `time_public_code.py`, which runs
inside X-PSI's own environment, imports X-PSI there, and builds the example's objects with lines copied from the
example's `main.py` (its docstring says so). X-PSI's own three tests of its two shipped likelihood values and of its
likelihood step pass on this laptop, and at its shipped settings it gives its shipped values here, −32552.616905 and
−47881.278177 (`CLAIMS.md`, H1 and H2). `compare_with_public_code.py` sets our programs beside what it saved; what that
printed on 23 November is in the folder `measured-2026-11-23/`. The setup page (`README.md` in this folder, section 6)
says how to repeat it.

### 1.3 Its time on this machine: measured on 23 November 2026

One likelihood of the quick example, a different star at every evaluation (12 stars, each number of the example's star
moved by up to 1 part in 1000), on one processor core of the same laptop (Intel Core Ultra 7 258V), one thread, with the
laptop on mains power (`CLAIMS.md`, H4 to H7):

| program and setting | time for one likelihood | X-PSI's time / ours | largest difference in ln L from ours converged, over the 12 stars |
|---|---|---|---|
| X-PSI 3.5.0b1 at its shipped settings (32 cells, 64 rotation steps, 512 rays, 64 energies) | fastest 49.1 ms, median 55.4 ms | 1 | 0.455 |
| ours compiled, 10 rings per panel, 16 points, 24 energies | 0.423 ms | 116 (against its fastest) to 131 (against its median) | 0.011 |
| ours compiled, 8 rings per panel, 12 points, 16 energies | 0.302 ms | 163 to 184 | 0.072 |
| ours in Python, 10 rings per panel, 16 points, 24 energies | fastest 2.94 ms, median 3.31 ms | 17 (fastest against fastest, and median against median) | 0.011 |
| ours in Python, 8 rings per panel, 12 points, 16 energies | fastest 2.41 ms, median 2.61 ms | 20 (fastest), 21 (median) | 0.072 |

So, compiled, our program works out the quick example's likelihood more than 100 times faster than X-PSI at X-PSI's
shipped settings, on one core of the same laptop: 116 to 131 times on mains power, with ours 0.011 from its converged value
and X-PSI 0.455 from it at the furthest of the 12 stars. In Python, 17 times. The compiled time is the mean over 300
different stars, the fastest of 7 rounds; in Python ours was timed as X-PSI was, one evaluation per star, the same 12 stars
five times over.

The same comparison on battery power, with X-PSI timed about ten minutes before ours: X-PSI fastest 73.5 ms, median
84.3 ms; ours compiled 0.652 ms at (10, 16, 24), 113 to 129 times, and 0.415 ms at (8, 12, 16), 177 to 203 times; in
Python fastest 4.17 ms and median 5.19 ms at (10, 16, 24), 18 times against its fastest and 16 against its median.

The last column is measured from ours converged (48 rings per panel, 32 points, 48 energies): three times finer settings
change ours by 0.00008 in ln L at this star and no bin by more than 4 parts in 10 million (H11), and as X-PSI's settings are
refined its counts move toward ours, per phase bin with its mesh and per channel with its energies (H8). X-PSI's own check
of this example allows 0.3255 at the example's star, where its value is 0.093 from ours converged.

At finer settings X-PSI takes longer: on mains power, 0.54 s at (64 cells, 128 rotation steps, 512 rays, 128 energies),
4.2 s at (128, 256, 512, 128), and 22.5 s at (192, 384, 512, 192) from a single evaluation (H10). Its score at the example's
star does not approach ours steadily as its settings are refined (+0.093 at its shipped settings, then −0.064, −0.080,
+0.029, +0.013 and +0.002 at finer ones; H10), because its errors in energy and in its mesh partly cancel. So no single
ratio "at equal accuracy" is stated: each ratio above is at X-PSI's shipped settings, with the accuracy of each program
beside it.

The times are for the quick example only. At the second shipped check the stopwatch's nearby stars fall half below X-PSI's
lower limit of the radius and are refused at once, so no time is stated for that check (`COVERAGE.md`, section 6, item 13).

What was said before the measurement, kept as it was said. X-PSI's own two statements of time are quoted below exactly as
its documentation prints them. Neither is a measurement on this machine; the measurement is above.

1. `docs/source/overview.rst`: "In general it is necessary to run X-PSI posterior sampling/integration on distributed
   memory architectures (e.g., a cluster or a supercomputer) because likelihood evaluation times are slow, of
   :math:`\mathcal{O}(1)\,s`." (The last phrase is written for the documentation's typesetter and renders as "of O(1) s".)
2. `docs/source/Modeling.ipynb`, the printed output of the cell that runs
   `likelihood.check(None, [-2.6848565597e+04], 1.0e-6, physical_points=[p])` between `t = time.time()` and
   `print('time = %.3f s' % (time.time() - t))`: "time = 0.662 s". The setting of that cell, as the notebook gives it: two
   hot regions (`hot = HotRegions((primary, secondary))`), each with `num_leaves=100` (its rotation steps), `num_rays=200`
   and `sqrt_num_cells=32`; the likelihood with `num_energies=128` and `threads=1`; the phase interpolant set to `'Akima'`
   in the cell before; a star of eleven parameters that is the notebook's own, not the quick example's; the notebook's own
   machine.

The project's plan also holds an estimate of the public code's time for the quick example, made before the measurement
from the public code's source under assumptions stated with it: 0.1 to 0.66 s, which "would give" a ratio of "30 to 200 at
the Python's 3.2 ms and 250 to 1,650 at the compiled 0.4 ms". The plan labels it, and a second estimate beside it, as
estimates: "Both are estimates. The measurement ends the question." Measured on 23 November 2026: 49.1 to 55.4 ms on mains
power, and ratios of 17 in Python and 116 to 131 compiled, all below the estimate's ranges. A later estimate of 0.1 to 0.2 s
was too high as well (`COVERAGE.md`, section 6, item 12).

## 2. The graphics-card code of Zhou and Huang

Its `README.md` names it "gpu_ppm: GPU-Accelerated Pulsar / Neutron Star Pulse Profile Modeling" and describes it as
"High-performance C++20 / CUDA code for generating phase-resolved pulse profiles and energy-phase counts for rotating
neutron star surface emission with relativistic light bending, time delays, Doppler / aberration effects, oblate star
shape, atmosphere beaming (NSX tables), and instrument response folding. Includes CPU reference implementations for
validation and comparison." Its paper is arXiv:2510.07764. Its licence file names Tianzhe Zhou and Chun Huang
(`LICENCE_NOTES.md`, section 4).

Its figures, read from the legends of the paper's figures (arXiv:2510.07764, version 1; its abstract says 2 to 5
thousandths of a second): 1.8 to 5.3 thousandths of a second per profile
at its standard setting and 6.9 to 14.2 thousandths of a second at its high setting, on one RTX 4080 graphics card; with a
hot-gas (atmosphere) spectrum, interstellar absorption and the instrument response, and without the likelihood.

It is the fastest public solution of this problem known to us.

The comparison with ours is not like for like, for these reasons:

- ours runs on one core of a laptop processor; theirs on a graphics card;
- ours includes the likelihood; theirs does not;
- ours uses a blackbody spectrum, no absorption, and the public example's instrument carried as one constant effective
  area folded with the energy integrals into one fixed matrix; theirs a hot-gas spectrum, absorption and the instrument
  response;
- the two were not run on the same machine, nor on the same problem.

Its source is in the reference material on this machine, but it was not run, and it is not needed here; only its
`README.md` and `LICENSE` were read.

## 3. Our own numbers on this page

One whole working-out of the public quick example's likelihood (the star's bin counts, the energy integrals, the
likelihood) at 10 rings per panel, 16 points and 24 energies, on one core of an Intel Core Ultra 7 258V, one process, one
thread, fastest of repeated runs, on 22 November 2026: 3.2 thousandths of a second in Python with NumPy, and 0.4 thousandths
of a second in the compiled copy. A repeat on 23 November with the laptop on battery power gave 5.1 and 1.0 thousandths of a
second, and on mains power again that day 2.9 and 0.40. Each
has its row in `CLAIMS.md` (F1, F2, F8) with the command that measured it and its date. Our likelihood values, and what they
do and do not show, are in `COVERAGE.md`, sections 3.1 and 5.

Beside X-PSI on 23 November 2026, the comparison timed ours as X-PSI was timed, at a different star every evaluation: on
mains power 0.423 ms compiled and fastest 2.94 ms in Python at (10, 16, 24), 0.302 ms and 2.41 ms at (8, 12, 16); on
battery power 0.652 and 0.415 ms compiled, 4.17 and 3.31 ms in Python (`CLAIMS.md`, H5 and H6; section 1.3 above).

## 4. What is not stated here

No ratio of our time to X-PSI's other than those of section 1.3, which are for the quick example, at X-PSI's shipped
settings, on one core of the same laptop, with the power state and the accuracy of each program beside each ratio. No
single ratio "at equal accuracy" (section 1.3 says why). No ratio for other problems (several spots, a tabulated
atmosphere: `COVERAGE.md`, section 4) or for other machines. No ratio to the graphics-card code's time:
not the same machine, not the same problem; nothing here says that ours is faster than it. No part of this work was asked
for, reviewed or endorsed by the authors of either code, by NASA or by the NICER team.
