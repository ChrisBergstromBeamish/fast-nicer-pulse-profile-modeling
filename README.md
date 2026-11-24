# Fast pulse profile modeling for NASA's NICER neutron-star observations

This repository holds a fast calculation of the X-ray pulse profile of a spinning neutron star with a hot spot, and of the
likelihood of observed photon counts given that profile. That calculation sits at the centre of pulse profile modeling of
NICER observations, where Bayesian inference over very many such likelihoods constrains the star's mass and radius, so its
speed sets the cost of the whole inference. The calculation is written in Python with NumPy so that every step can be read
and checked, and written out again in C for speed; both are compared, bin by bin and in the likelihood, with X-PSI, the
public code for this problem, on X-PSI's own shipped example. The only data used so far is that synthetic example.

NASA, the NICER team and the authors of the public codes have not asked for, reviewed or endorsed this work. NICER is named
here only as the instrument whose kind of data this calculation serves.

## The three folders

1. `1-setup-and-documentation/`: how to set the repository up (`README.md`), what the programs compute and from what
   (`DERIVATION.md`), every stated number with the command that reproduces it (`CLAIMS.md`), what is covered and what is
   not (`COVERAGE.md`), the public codes (`PUBLIC_CODES.md`), the licences read (`LICENCE_NOTES.md`), the roadmap, and
   what the comparison with the public code X-PSI printed on 23 November 2026 (`measured-2026-11-23/`).
2. `2-python-solution/`: the solution in Python with NumPy. This is the reference: it is the one that can be read line by
   line against the derivation. It also holds a second, slow, independent program that the fast one is checked against, and
   the one command that repeats every check (`check.py`).
3. `3-compiled-solution/`: the same steps written out again in C, built in two forms by one `make`, and checked against the
   Python solution bin by bin by one command (`check.py`).

## What is measured so far

One processor core of a laptop (Intel Core Ultra 7 258V), unless a line says otherwise; every number below has a row in
`1-setup-and-documentation/CLAIMS.md` with the command that reproduces it.

- The likelihood of the public code X-PSI's shipped quick example: -32552.710 against its shipped -32552.617, where its
  own check accepts a difference of 0.3255.
- One working-out of that likelihood: 3.2 thousandths of a second in Python, 0.4 thousandths of a second compiled, as
  measured on 22 November 2026. A repeat on 23 November with the laptop on battery power gave 5.1 and 1.0 thousandths of
  a second with the same numbers out; on mains power again that day, 2.9 and 0.40. All are in the claims table.
- Many stars at once, on all eight processors of the same laptop (sprint 2): on eight threads (lines of work of one
  program, each on a processor core) the compiled program worked out 2,096 stars per second at the working setting
  (10 rings per panel, 16 points, 24 energies; median 2,086), 5.2 times its 406 per second on one thread, on batches of
  4,096 different stars drawn from the public example's prior, fastest of five repeats, on mains power with nothing else
  running on 24 November 2026 (`CLAIMS.md`, row I1). In the checks, each of 108 stars came out of batches on 1 to 8 threads the same, bit for
  bit, as when it is worked out alone (row J1). No ratio against X-PSI is stated for several cores: X-PSI's time on several threads is
  not measured.
- The public code X-PSI itself (version 3.5.0b1), installed in its own folder on the same laptop and run there on
  23 November 2026: its own tests pass, and at its shipped settings it gives its shipped values. Compiled, our program
  works out the quick example's likelihood more than 100 times faster than X-PSI at X-PSI's shipped settings, on one core
  of the same laptop: 116 to 131 times on mains power (X-PSI 49 to 55 thousandths of a second, ours 0.42), 113 to 129 times
  on battery. In Python, 17 times on mains power. Beside the ratios, the accuracy: over twelve stars near the example's,
  ours at that setting is at most 0.011 from our converged value and X-PSI at its shipped settings at most 0.455; as
  X-PSI's mesh and its energies are refined its counts move toward ours (its rotation steps and its rays alone change
  little), and ours do not move when ours are refined. X-PSI's score does not approach ours steadily as its settings
  are refined, so no single ratio "at equal accuracy" is stated (`1-setup-and-documentation/PUBLIC_CODES.md`,
  section 1.3).
- The public code's second shipped check (a different star): ours gives +6.5 against a tolerance of 0.479, so ours fails
  that check as shipped. The cause is measured: the shipped number is X-PSI at coarse settings for that star. With its
  settings refined (128 cells, 256 rotation steps, 512 rays, 384 energies) X-PSI itself gives -47874.997, 6.28 above its
  shipped number and 0.24 from ours, -47874.760 (`COVERAGE.md`, section 5.3).

Agreement with the public code is not a result about the physics, and never a proof of it: `DERIVATION.md` says from which
starting points the programs compute light's path, and why this project cannot tell those starting points apart from the
public code's. The result of this project is a method of calculation.

## Licence

GNU General Public License, version 3: the file `LICENSE`. `1-setup-and-documentation/LICENCE_NOTES.md` says which public sources were read and under which licences they stand.

## What is not in this repository

The public code X-PSI, its example data and the published reference curves are not copied here. The setup page says where
they come from and where to put them. The checks do not need X-PSI installed; the comparison of times and bins with it does,
and the scripts that install and run it for that comparison are in `2-python-solution/tools/` (setup page, section 6).
