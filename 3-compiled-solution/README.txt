The compiled copy of the fast program

nicer_fast.c      the Python program written out again in C: the star's bin counts
                  (2-python-solution/code/fast_kernel.py, spot_counts, with the one fit of light's path) and the public
                  likelihood's definition (2-python-solution/code/public_likelihood.py). Nothing of the public code X-PSI
                  is in it. The Python programs are the reference.
tables.bin        what is made once and never changes between workings-out. Not kept in the repository; written by
                      python -I -B -X utf8 ../2-python-solution/tools/export_for_c.py <the public example's Data folder> tables.bin
time_compiled.sh  times the program at several settings.

What it is compiled into: one native program for the processor it is built on (here: 64-bit Intel/AMD, Linux).
What it needs: a C compiler (C99) and the standard C and math libraries. Nothing else: no Python at run time, no other
library. It was built and timed on Ubuntu 24.04 (under Windows' Linux layer, WSL) with gcc 13.3.

Build:
    gcc -O3 -march=native -ffast-math -Wall -Wextra -o nicer_fast nicer_fast.c -lm
Check against the Python, bin by bin, on six stars at three settings (run on the Windows side; it calls the Linux side):
    python -I -B -X utf8 ../2-python-solution/tools/compare_c_with_python.py <the public example's Data folder>

Arguments: nicer_fast <tables> <rings per panel> <points> <energies> <repeats> mass radius distance cos_i phase_shift colat rho log10T [dump file]
    repeats > 0: time that many workings-out of the same star; repeats < 0: a different star at each one; 0: the value only.
It prints "refused: code -1" for a star outside the range it covers (L/R between 0.15 and 0.55, moderate flattening),
-2 where light from behind the star would arrive by a second path, -3 if an end direction of a ring was not found.

Measured 22 Nov 2026, Intel Core Ultra 7 258V, one thread: 0.4 ms per working-out at 10 rings per panel, 16 points,
24 energies (ln L = -32552.719298 for the public quick example; the Python gives the same to 2e-9 and takes 3.2 ms).
