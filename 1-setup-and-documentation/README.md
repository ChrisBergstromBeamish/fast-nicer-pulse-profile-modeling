# Setup and documentation

Seed page (23 November 2026). Sprint 1 completes it.

## What the programs need

- The Python solution: Python 3.12 and NumPy. Nothing else.
- The compiled solution: a C compiler (C99) and the standard C and math libraries. Built so far with gcc 13.3 on
  Ubuntu 24.04. It reads one table file that the Python solution writes once.

## Public material the checks read (not copied into this repository)

- X-PSI, commit bb5f461e605b22c6251b155b27a914a558757f93: the data files of its quick example
  (`examples/examples_fast/Data`) and the expected values in its tests. Only read; X-PSI is not needed to run our
  programs.
- The reference curves and tables of Bogdanov et al. 2019 (Astrophysical Journal Letters 887, L26), supplementary
  data `apjlab5968`.

## The quick checks, as they run today

    python -I -B -X utf8 2-python-solution/tools/run_all_checks.py <X-PSI's examples/examples_fast/Data folder> <the apjlab5968 folder>

## Still to be written (sprint 1)

- The derivation, step by step: what is assumed, what is derived, what is taken from the public example as its input.
- For every stated number: the command that reproduces it and the value to expect.
- What the programs cover and where they refuse.
