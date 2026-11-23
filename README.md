# Fast pulse profiles

A fast calculation of the X-ray pulse profile of a spinning neutron star with a hot spot, and of the likelihood of
observed counts given that profile, as used for NICER data.

Status: seed copy (23 November 2026). The programs are as they were first written and tested. Sprint 1 makes the
repository stand on its own: one way to set it up, one command that repeats every check, and documentation that shows
the work step by step.

## The three folders

1. `1-setup-and-documentation/`: how to set it up, what each program computes and from what, and how to check every
   stated number yourself.
2. `2-python-solution/`: the solution in Python with NumPy. This is the reference: it is the one that can be read line
   by line against the derivation. It also holds a second, slow, independent program that the fast one is checked
   against.
3. `3-compiled-solution/`: the same steps written out again in C and compiled. It is checked against the Python
   solution bin by bin.

## What is measured so far (one processor core of a laptop, Intel Core Ultra 7 258V)

- The likelihood of the public code X-PSI's shipped quick example: -32552.710 against its shipped -32552.617, where its
  own check accepts a difference of 0.3255.
- One working-out of that likelihood: 3.2 thousandths of a second in Python, 0.4 thousandths of a second compiled.
- The public code's time on the same laptop: not measured yet.

## What is not in this repository

The public code X-PSI, its example data and the published reference curves are not copied here. The setup page says
where they come from.
