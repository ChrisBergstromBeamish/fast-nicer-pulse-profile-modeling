# The compiled solution

The fast program and the public likelihood's definition, written out again in C, in one file. The Python solution is the
reference; this folder's program is checked against it bin by bin (`check.py`). Nothing of the public code X-PSI is in it.

| file | what it is |
|---|---|
| `nicer_fast.c` | the Python program written out again in C: the star's bin counts (`2-python-solution/code/fast_kernel.py`, `spot_counts`, with the approximation of the path integrals) and the public likelihood's definition (`2-python-solution/code/public_likelihood.py`) |
| `Makefile` | `make` writes the table file and builds the program in two forms; `make depends` reports what was built; `make clean` |
| `check.py` | one command that compares both forms with the Python solution (run from the repository's top, on the side that has Python with NumPy) |
| `time_compiled.sh` | times the program at several settings (prints the processor and compiler first) |
| `tables.bin` | made on each machine by `make`, never kept: everything that never changes between workings-out, written by `2-python-solution/tools/export_for_c.py` from our own Python programs and the public example's data |
| `nicer_fast`, `nicer_fast_strict` | the two built forms, made on each machine, never kept |

## What it is compiled into, and what it needs to run

Both forms are native programs for the processor they are built on: `-march=native` lets the compiler use every instruction
of that processor, so a built program runs on that machine and on others with the same instruction extensions, and is
rebuilt (one `make`) for any other. Found with `make depends`, which runs `file` and `ldd` on the two programs, on the
machine the work was done on (Ubuntu 24.04 under WSL, gcc 13.3, a 64-bit Intel processor):

- `file`: each form is an "ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter
  /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=..., for GNU/Linux 3.2.0, ..., not stripped": a 64-bit Linux program for
  Intel/AMD processors that loads shared libraries at start. For `nicer_fast` the line also says "with debug_info", for
  `nicer_fast_strict` it does not; neither is built with `-g`, and the difference was not traced further.
- `ldd nicer_fast` (the compiler's fast arithmetic): `libm.so.6` (the C math library), `libmvec.so.1` (the C library's
  vector math library: with `-ffast-math` the compiler calls its vectorised sin, cos, exp and log), `libc.so.6` (the C
  library), `/lib64/ld-linux-x86-64.so.2` (the dynamic loader) and `linux-vdso.so.1` (the kernel's own, not a file).
- `ldd nicer_fast_strict` (strict arithmetic): the same without `libmvec.so.1`.

So at run time each form needs: a 64-bit Linux with the GNU C library (which provides all of the above; glibc 2.39 was
used, as `ldd --version` prints it; the compiler was gcc 13.3.0, as `gcc --version` prints it), the table file `tables.bin`
next to it, and nothing else: no Python, no other library. Python with NumPy is needed once, when `make` writes
`tables.bin`.

The two forms, built by `make` from the same source:

| form | compiler command | arithmetic |
|---|---|---|
| `nicer_fast` | `gcc -O3 -march=native -ffast-math -Wall -Wextra -o nicer_fast nicer_fast.c -lm` | the compiler may reorder floating-point operations, assume no infinities and no not-a-numbers, and fuse multiply-add; the form timed so far |
| `nicer_fast_strict` | `gcc -O3 -march=native -ffp-contract=off -Wall -Wextra -o nicer_fast_strict nicer_fast.c -lm` | every operation as written, in order, rounded as written: no `-ffast-math` and no fused multiply-add |

Both forms give the same numbers as the Python solution to the limits in `check.py` (the fast form was also compared with
the Python solution with the limits of `1-setup-and-documentation/CLAIMS.md`; the near-axis star of the six-star check is
the one case where the two differ by more than 1e-10 of the mean bin, in both forms alike).

## Running it by hand

    nicer_fast <tables> <rings per panel> <points> <energies> <repeats> mass radius distance cos_i phase_shift colat rho log10T [dump file]

`repeats > 0`: time that many workings-out of the same star; `repeats < 0`: a different star at each one (as a sampler asks);
`0`: the value only. The eight numbers are the public example's parameters in its own names and units (solar masses, km,
kpc, the cosine of the inclination, cycles, radians, radians, log10 of kelvin). With a dump file the bin counts and the
model's count rates are written to it as raw doubles (read by `check.py`'s scripts). A second mode,
`nicer_fast <tables> channels <file>`, runs the likelihood step alone on channels read from a file (used by the hard-channel
check).

It prints `refused: code -1` for a star outside the range it covers (L/R between 0.15 and 0.55, the flattening
`|o2|` below 0.5 and the surface's tilt below the approximation's bound), `-2` where light from behind the star could reach
the observer by a second path, `-3` if an end direction of a ring was not found; `1-setup-and-documentation/COVERAGE.md`
has the full list. `ln L = -inf` is the public code's own early exit (the model far above the data in some channel).

Times measured on this machine are in `1-setup-and-documentation/CLAIMS.md`, each with the command that measured it.
