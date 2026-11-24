# The compiled solution

The fast program and the public likelihood's definition, written out again in C, in one file. The Python solution is the
reference; this folder's program is checked against it bin by bin (`check.py`). Nothing of the public code X-PSI is in it.

Three terms of the trade used below: a thread is one line of work of a program, each on a processor core, so that several
run at the same moment; OpenMP is a feature of the compiler that shares the turns of a loop among threads; a race is two
threads touching the same memory with nothing to put one before the other.

**To build**, on the Linux side, in this folder: `make` (it writes `tables.bin` with a Python that has NumPy, then builds
the four forms below; the setup page, `1-setup-and-documentation/README.md`, section 5, says what it needs and how it
finds that Python). `make race` builds and runs the race detector's form, `make depends` reports what was built, and
`make clean` removes it all. Then `check.py` runs every check, from the repository's top.

| file | what it is |
|---|---|
| `nicer_fast.c` | the Python program written out again in C: the star's bin counts (`2-python-solution/code/fast_kernel.py`, `spot_counts`, with the approximation of the path integrals) and the public likelihood's definition (`2-python-solution/code/public_likelihood.py`) |
| `Makefile` | `make` writes the table file and builds the program in two forms, the fast form without OpenMP, and the library; `make race` runs the race detector; `make depends` reports what was built; `make clean` |
| `check.py` | one command that compares both forms with the Python solution, and the batch mode with each star alone (run from the repository's top, on the side that has Python with NumPy) |
| `time_compiled.sh` | times the program at several settings on one thread (prints the processor and compiler first) |
| `time_batch.sh` | times the batch mode on 1, 2, 4, 6 and 8 threads (prints the processor, the compiler and the power state) |
| `nicer_fast_lib.py` | calls the library from Python with only Python's standard library (`ctypes`) |
| `draw_stars.py` | writes stars drawn from the public example's prior, one per line, for the batch mode (standard library only) |
| `each_alone.sh`, `batch_check.sh` | run on the Linux side by `check.py`: each star of a file alone, then the same stars in batches |
| `tables.bin` | made on each machine by `make`, never kept: everything that never changes between workings-out, written by `2-python-solution/tools/export_for_c.py` from our own Python programs and the public example's data |
| `nicer_fast`, `nicer_fast_strict`, `nicer_fast_no_openmp`, `libnicer_fast.so` | the built forms, made on each machine, never kept (also `nicer_fast_race`, made by `make race`) |

## What it is compiled into, and what it needs to run

All built forms are native code for the processor they are built on: `-march=native` lets the compiler use every instruction
of that processor, so a built program runs on that machine and on others with the same instruction extensions, and is
rebuilt (one `make`) for any other. Found with `make depends`, which runs `file` and `ldd` on the four built forms, on the
machine the work was done on (Ubuntu 24.04 under WSL, gcc 13.3, a 64-bit Intel processor):

- `file`: each of the three programs (`nicer_fast`, `nicer_fast_strict`, `nicer_fast_no_openmp`) is an "ELF 64-bit LSB pie executable, x86-64, version 1 (SYSV), dynamically linked, interpreter
  /lib64/ld-linux-x86-64.so.2, BuildID[sha1]=..., for GNU/Linux 3.2.0, ..., not stripped": a 64-bit Linux program for
  Intel/AMD processors that loads shared libraries at start. For `nicer_fast` the line also says "with debug_info", for
  `nicer_fast_strict` it does not; neither is built with `-g`, and the difference was not traced further.
- `ldd nicer_fast` (the compiler's fast arithmetic): `libm.so.6` (the C math library), `libmvec.so.1` (the C library's
  vector math library: with `-ffast-math` the compiler calls its vectorised sin, cos, exp and log), `libc.so.6` (the C
  library), `/lib64/ld-linux-x86-64.so.2` (the dynamic loader) and `linux-vdso.so.1` (the kernel's own, not a file).
- `ldd nicer_fast` also lists `libgomp.so.1`: gcc's OpenMP library, which runs the batch mode's threads. OpenMP is a
  feature of the compiler that shares the turns of a loop among threads (lines of work of one program, each on a core);
  `-fopenmp` turns it on, and adds this one library.
- `ldd nicer_fast_strict` (strict arithmetic): the same without `libmvec.so.1`.
- `nicer_fast_no_openmp`: the same `file` line as `nicer_fast`; `ldd` lists what `nicer_fast` lists without `libgomp.so.1`.
- `libnicer_fast.so`: `file` says "ELF 64-bit LSB shared object, x86-64, version 1 (SYSV), dynamically linked,
  BuildID[sha1]=..., not stripped"; `ldd` lists the same libraries as `nicer_fast`.

So at run time each form needs: a 64-bit Linux with the GNU C library (which provides all of the above but `libgomp.so.1`;
glibc 2.39 was used, as `ldd --version` prints it; the compiler was gcc 13.3.0, as `gcc --version` prints it), gcc's
`libgomp.so.1` for every form but `nicer_fast_no_openmp`, the table file `tables.bin` next to it, and nothing else: no
Python, no other library. Python with NumPy is needed once, when `make` writes `tables.bin`; the library's caller
`nicer_fast_lib.py` needs a Python with its standard library only.

The forms, built by `make` from the same source:

| form | compiler command | arithmetic |
|---|---|---|
| `nicer_fast` | `gcc -O3 -march=native -ffast-math -Wall -Wextra -fopenmp -o nicer_fast nicer_fast.c -lm` | the compiler may reorder floating-point operations, assume no infinities and no not-a-numbers, and fuse multiply-add; the form timed so far |
| `nicer_fast_strict` | `gcc -O3 -march=native -ffp-contract=off -Wall -Wextra -fopenmp -o nicer_fast_strict nicer_fast.c -lm` | every operation as written, in order, rounded as written: no `-ffast-math` and no fused multiply-add |
| `nicer_fast_no_openmp` | `gcc -O3 -march=native -ffast-math -Wall -Wextra -o nicer_fast_no_openmp nicer_fast.c -lm` | as `nicer_fast`; the batch mode runs on one thread whatever number is asked for |
| `libnicer_fast.so` | `gcc -O3 -march=native -ffast-math -Wall -Wextra -fopenmp -fPIC -shared -DNICER_LIBRARY -o libnicer_fast.so nicer_fast.c -lm` | as `nicer_fast` |

`-fopenmp` was added in sprint 2 and is not an arithmetic option: the arithmetic options of the two forms are those they had
before, and every star gives the same bits as before (below).

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
has the full list. `ln L = -inf` is the public code's own early exit (the model far above the data in some channel). After
the line `ln L = ...` it prints `every bit of ln L: ...`, the same number as a hexadecimal floating-point number (C's `%a`),
which loses no bit.

## Many stars at once: the batch mode

A sampler asks for many likelihoods at a time, one per star it wants to try. The batch mode works out many stars at once,
one star per thread at any moment; inside one star nothing runs in parallel, and each star is worked out by the same
function, in the same order, as alone.

    nicer_fast <tables> batch <rings per panel> <points> <energies> <threads> <stars file> [dump file]

The stars file holds one star per line, the eight numbers in the order above, separated by spaces (blank lines are skipped).
Numbers after the eighth on a line are not read. The batch mode uses the example's limits on the background
(`NICER_NO_UPPER_LIMIT` is read only by the single-star mode). For each star, in the file's order, one line is written: ln L with up to 17 significant digits (C's `%.17g`, which leaves out trailing zeros, so a line may show 16 or fewer; read as a double it gives back every bit
as a double), or `-inf` for the public code's early exit and for a refused star, then the code (0, or the refusal's code as
above; `-9` if a thread could not have its working memory). A refused star does not disturb its neighbours. With a dump
file, each star's bin counts and model rates are written to it, one star after another (zeros for a refused star). The time
the batch took, and stars per second, go to the error output (the stars only, not the start of the program or the loading
of the tables). For example, with 256 stars drawn from the public example's prior:

    python3 -I -B draw_stars.py 1 256 stars_example.txt
    ./nicer_fast tables.bin batch 10 16 24 8 stars_example.txt > batch_example.txt

(`make clean` removes both files; `.gitignore` lists `stars_*.txt` and `batch_*.txt`.)

The stars are handed to the threads one at a time as threads become free (OpenMP's `schedule(dynamic, 1)`), since stars
differ in how long they take. `nicer_fast_strict` has the same batch mode; `nicer_fast_no_openmp` runs it on one thread and
says so if more are asked for.

**From Python.** `libnicer_fast.so` is the same code as a shared library with three entry points: `nicer_load(tables)` (once
per process; a second call keeps the tables of the first; the settings too are the process's, set by the last
`nicer_choose`), `nicer_choose(rings per panel, points, energies)` and `nicer_batch(N, stars, threads, lnL, codes, dump or
NULL)`. `nicer_fast_lib.py` calls them with Python's standard library alone:

    import nicer_fast_lib
    lib = nicer_fast_lib.Batch("tables.bin", 10, 16, 24)
    lnl, codes = lib.work_out(stars, threads=8)        # stars: a list of stars, each eight numbers; -1e300 stands for -inf

Run as a command (`python3 -I -B nicer_fast_lib.py tables.bin 10 16 24 8 stars_example.txt`) it prints the batch mode's
lines; `check.py` compares the two. A missing or wrong table file ends the calling process with a message, as it ends the
program.

**The checks of the batch mode** (in `check.py`, after the comparisons with Python). Each form gives every star in a batch
exactly the bits it gets alone (ln L, the code, every bin count and every model rate): the sixty prior draws, the six stars at
three settings, and ten stars among which six are refused (code -1) and three take the early exit, at three settings; on 1, 2,
3, 4 and 8 threads, and again with the stars in a shuffled order. The same batch run twice gives the same bits. The library
gives the batch mode's lines. And `make race` builds `nicer_fast_race` with the compiler's race detector (`-fsanitize=thread`:
it watches a running program for two threads touching the same memory with nothing to put one before the other), runs 256
stars on 8 threads, and fails if it reports a race. It needs the detector's library `libtsan` (part of gcc on Ubuntu 24.04)
and runs the program under `setarch -R` (no address randomisation), without which this detector stops at once on recent
Linux kernels ("unexpected memory mapping"). The OpenMP library is not built for the detector, so the detector cannot see
the barrier at the end of the batch's loop, where OpenMP makes the main thread wait until every thread has written its
stars; it first reported the main thread's reading of the results as a race. `nicer_batch` therefore ends with a counter
that each thread raises after its writes and that the main thread reads after the loop (an atomic "release" and "acquire",
which the detector does see); it orders nothing that the barrier did not already order, and changes no number. Nothing of
the OpenMP library itself was reported, so `make race` leaves nothing out of the detector's watch (its option
`ignore_noninstrumented_modules` is not used). That the checks can fail was shown once, in a copy that is not kept: with one
working buffer shared by every thread, 49 of the sixty stars differed from their values alone on 8 threads (none on one
thread; a second run found 55: the count varies from run to run, as a race's does), and the race detector
stopped at the first shared write, in `spot_counts`. `make race` runs one batch in one process; a second batch in the same
process (as a Python sampler would ask) is not run under the detector.

Before the batch mode was added, both forms were built at the commit before it and their outputs kept, for the quick example
at (10, 16, 24), (24, 24, 32) and (8, 12, 16) and for the sixty, six and ten stars above at the same three settings: every
bin count of the dump and every bit of ln L. The forms built now give the same bytes alone and in every batch, on 1 to 8
threads (that comparison was made once, in a folder that is not kept; `check.py` repeats the part that compares the batch
with each star alone).

Times measured on this machine are in `1-setup-and-documentation/CLAIMS.md`, each with the command that measured it. The
batch mode's times on 1 to 8 threads (`time_batch.sh`) are its rows I1 to I7.

## What a working-out writes, and where

Two workings-out that run at the same moment on two threads (two lines of work of one program, each on a core) must never
write to the same memory. Every variable of `nicer_fast.c` that lives outside a function, or keeps its value between calls
(`static`), is in this table: set once before any working-out and only read afterwards, or written during a working-out.

| name | where | which kind | what was done |
|---|---|---|---|
| `arrs`, `narr` | file level | the tables: set once by `load_tables` (`nicer_load`) | kept; during a working-out only read (`get`, `getn`) |
| `NU`, `NZ`, `ND`, `NVIS`, `NPP`, `NN`, `NE`, `NB`, `NCH` | file level | the settings: set once by `choose` (`nicer_choose`); the `channels` mode sets `NCH`, `NB` once in `main` | kept; only read |
| `FIT_CT`, `FIT_CA`, `VUG`, `VT0`, `VT0X`, `VT0XX`, `GLG`, `GLW`, `CCOS`, `CSIN`, `CG`, `CQ0`, `CTOD`, `EN_E`, `EN_W`, `COUNTS`, `SUPPORT`, `FACT` | file level | pointers into the tables: set once by `choose`; `SUPPORT` also once by `main` (`NICER_NO_UPPER_LIMIT`), and the `channels` mode sets `COUNTS`, `SUPPORT`, `FACT` once | kept; only read |
| `FIT_U0`, `FIT_U1`, `ZM`, `FIT_XMIN`, `FIT_TILT`, `C_CM`, `XGM`, `XKPC`, `XKB`, `XBB`, `T_EXP`, `FREQ` | file level | numbers from the tables: set once by `choose` (`T_EXP` by the `channels` mode) | kept; only read |
| `GQ`, `WQ` | file level | the likelihood's quadratures: set once (`nicer_choose`; the `channels` mode) | kept; only read |
| `NQ`, `QLIM` | file level | constants (`const`) | never written |
| `PER_CHANNEL` | file level (a pointer) | the array it pointed to was written during a working-out (each channel's value, for the test of hard channels) | removed from file level: now an argument of `marginal_loglike`, owned by the caller (`NULL` in a working-out of a star) |
| `theta`, `wth`, `dlam` | `static` inside `spot_counts` | written during every working-out (the rings of the spot) | moved to the function's own stack |
| `S`, `star` | `static` inside `main` | written during every working-out of the single-star mode | now taken by `main` with `malloc`; in a batch each thread takes its own `S` and `star` |

No other function keeps anything between calls. In a batch the tables and the settings are shared and only read.
