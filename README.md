# Fast pulse profile modeling for NASA's NICER neutron-star observations

This repository holds a fast calculation of the X-ray pulse profile of a spinning neutron star with a hot spot, and of the
likelihood of observed photon counts given that profile. That calculation sits at the centre of pulse profile modeling of
NICER observations, where Bayesian inference over very many such likelihoods constrains the star's mass and radius, so its
speed sets the cost of the whole inference. The calculation is written in Python with NumPy so that every step can be read
and checked, and written out again in C for speed; both are compared, bin by bin and in the likelihood, with X-PSI, the
public code for this problem, on X-PSI's own shipped example. The only data used so far is that synthetic example.

Suggested GitHub description: Fast pulse profile modeling for NICER: the X-ray pulse profile of a hot spot on a spinning neutron star and its likelihood, in Python and in C, checked bin by bin and against X-PSI's shipped example.

Suggested GitHub topics: nicer, pulse-profile-modeling, neutron-star, x-ray-astronomy, hot-spot, likelihood, bayesian-inference, x-psi, light-bending, python, c, numpy

NASA, the NICER team and the authors of the public codes have not asked for, reviewed or endorsed this work. NICER is named
here only as the instrument whose kind of data this calculation serves.

## The three folders

1. `1-setup-and-documentation/`: how to set the repository up (`README.md`), what the programs compute and from what
   (`DERIVATION.md`), every stated number with the command that reproduces it (`CLAIMS.md`), what is covered and what is
   not (`COVERAGE.md`), the public codes (`PUBLIC_CODES.md`), the licences read (`LICENCE_NOTES.md`), the roadmap.
2. `2-python-solution/`: the solution in Python with NumPy. This is the reference: it is the one that can be read line by
   line against the derivation. It also holds a second, slow, independent program that the fast one is checked against, and
   the one command that repeats every check (`check.py`).
3. `3-compiled-solution/`: the same steps written out again in C, built in two forms by one `make`, and checked against the
   Python solution bin by bin by one command (`check.py`).

## What is measured so far

One processor core of a laptop (Intel Core Ultra 7 258V); every number below has a row in
`1-setup-and-documentation/CLAIMS.md` with the command that reproduces it.

- The likelihood of the public code X-PSI's shipped quick example: -32552.710 against its shipped -32552.617, where its
  own check accepts a difference of 0.3255.
- One working-out of that likelihood: 3.2 thousandths of a second in Python, 0.4 thousandths of a second compiled, as
  measured on 22 November 2026. A repeat on 23 November on the same machine gave 5.1 and 1.0 thousandths of a second with
  the same numbers out; both are in the claims table.
- The public code's time on the same laptop: not measured. So no speed ratio against it is stated anywhere here.
- The public code's second shipped check (a different star): ours gives +6.5 against a tolerance of 0.479, so ours fails
  that check; which of the two numbers is nearer the exact value of the model is not known (`COVERAGE.md`).

Agreement with the public code is not a result about the physics, and never a proof of it: `DERIVATION.md` says from which
starting points the programs compute light's path, and why this project cannot tell those starting points apart from the
public code's. The result of this project is a method of calculation.

## What is not in this repository

The public code X-PSI, its example data and the published reference curves are not copied here. The setup page says where
they come from and where to put them.
