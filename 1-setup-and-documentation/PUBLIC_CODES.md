# The public codes this work is set beside

Two public codes solve the same kind of problem as the programs in this repository: X-PSI, the public code this work is
compared with, and the graphics-card code of Zhou and Huang. This page says what each is, what of it was used here and what
was not, what each says of its own speed, and why neither comparison with our programs is like for like yet. NASA, the
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

Not used: its code. X-PSI has not been installed, built or run on this machine; the decision to install it is the
author's. Nothing of it is imported by anything here, and no file of it
is copied here (`LICENCE_NOTES.md`). What the agreement with its shipped values does and does not show is in
`COVERAGE.md`, section 5.

### 1.3 Its time on this machine: not measured

Because the public code has not run here, its time for the quick example on this machine is not measured, and no speed
ratio of our programs against it may be stated as a fact. Its own two statements of time are quoted below exactly as its
documentation prints them. Neither is a measurement on this machine.

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

The project's plan also holds an estimate of the public code's time for the quick example, made from the public code's
source under assumptions stated with it: 0.1 to 0.66 s, which "would give" a ratio of "30 to 200 at the Python's 3.2 ms and
250 to 1,650 at the compiled 0.4 ms". The plan labels it, and a second estimate beside it, as estimates: "Both are
estimates. The measurement ends the question." Nothing on this page turns either into a fact.

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

## 4. What is not stated here

No ratio of our time to X-PSI's: X-PSI's time on this machine is not measured. No ratio to the graphics-card code's time:
not the same machine, not the same problem. No part of this work was asked for, reviewed or endorsed by the authors of
either code, by NASA or by the NICER team.
