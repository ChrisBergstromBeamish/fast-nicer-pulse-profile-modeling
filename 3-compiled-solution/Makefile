# The compiled solution: builds the table file and the program in its forms. Run on the Linux side, in this folder:
#     make            writes tables.bin (once), builds nicer_fast (the compiler's fast arithmetic, the form timed so far),
#                     nicer_fast_strict (strict arithmetic: no -ffast-math, and no fused multiply-add either), both with
#                     OpenMP for the batch mode; nicer_fast_no_openmp (the fast form built without OpenMP: the batch mode on
#                     one thread); and libnicer_fast.so (the fast form as a shared library, called by nicer_fast_lib.py)
#     make race       builds nicer_fast_race with the compiler's race detector and runs a batch of 256 stars on 8 threads
#     make depends    says what each built program is and which libraries it loads (file, ldd)
#     make clean      removes what make made
# Then, from the repository's top, on the side that has Python with NumPy:  python -I -B -X utf8 3-compiled-solution/check.py
# tables.bin is written by the Python solution (../2-python-solution/tools/export_for_c.py) and needs a Python with NumPy. The
# first of python3, python and python.exe that can import NumPy is used; on the Linux side of a Windows machine whose Linux
# Python has no NumPy, python.exe is the Windows Python, reached through WSL. Or name one:  make PYTHON=/path/to/python
# draw_stars.py (make race) needs only Python's standard library: PYTHON3 (python3 unless named) runs it.
# The public example's data folder comes from ../settings.cfg (1-setup-and-documentation/README.md).
CC = gcc
CFLAGS_FAST = -O3 -march=native -ffast-math -Wall -Wextra
CFLAGS_STRICT = -O3 -march=native -ffp-contract=off -Wall -Wextra
OPENMP = -fopenmp
PYTHON ?= $(shell for p in python3 python python.exe; do if $$p -I -c "import numpy" >/dev/null 2>&1; then echo $$p; break; fi; done)
PYTHON3 ?= python3

all: tables.bin nicer_fast nicer_fast_strict nicer_fast_no_openmp libnicer_fast.so

nicer_fast: nicer_fast.c
	$(CC) $(CFLAGS_FAST) $(OPENMP) -o $@ $< -lm

nicer_fast_strict: nicer_fast.c
	$(CC) $(CFLAGS_STRICT) $(OPENMP) -o $@ $< -lm

nicer_fast_no_openmp: nicer_fast.c
	$(CC) $(CFLAGS_FAST) -o $@ $< -lm

libnicer_fast.so: nicer_fast.c
	$(CC) $(CFLAGS_FAST) $(OPENMP) -fPIC -shared -DNICER_LIBRARY -o $@ $< -lm

tables.bin: ../2-python-solution/tools/export_for_c.py ../2-python-solution/code/fast_kernel.py ../2-python-solution/code/quick_example.py ../2-python-solution/code/public_likelihood.py ../2-python-solution/code/consts.py
	@test -n "$(PYTHON)" || { echo "no Python with NumPy was found (tried python3, python, python.exe): run  make PYTHON=<a python with NumPy>" >&2; exit 1; }
	$(PYTHON) -I -B -X utf8 ../2-python-solution/tools/export_for_c.py $@

# The race detector (-fsanitize=thread) watches the running program for two threads touching the same memory without order.
# It needs the compiler's detector library (libtsan, part of gcc on Ubuntu). The OpenMP library itself is not built for the
# detector; nothing of it was reported, so no report is left out (README.md). setarch -R runs the program without address
# randomisation, which this detector needs on recent Linux kernels (without it, it stops with "unexpected memory mapping").
nicer_fast_race: nicer_fast.c
	$(CC) $(CFLAGS_FAST) $(OPENMP) -g -fsanitize=thread -o $@ $< -lm

stars_race.txt: draw_stars.py
	$(PYTHON3) -I -B draw_stars.py 2 256 $@

race: tables.bin nicer_fast_race stars_race.txt
	TSAN_OPTIONS="halt_on_error=1 exitcode=66" setarch $$(uname -m) -R ./nicer_fast_race tables.bin batch 10 16 24 8 stars_race.txt > race_out.txt
	@echo "the race detector found no race: 256 stars on 8 threads"

depends: nicer_fast nicer_fast_strict nicer_fast_no_openmp libnicer_fast.so
	file nicer_fast nicer_fast_strict nicer_fast_no_openmp libnicer_fast.so
	ldd nicer_fast
	ldd nicer_fast_strict
	ldd nicer_fast_no_openmp
	ldd libnicer_fast.so

clean:
	rm -f nicer_fast nicer_fast_strict nicer_fast_no_openmp libnicer_fast.so nicer_fast_race tables.bin hard_channels*.bin hard_channels*.bin.out dump_*.bin stars_*.txt race_out.txt batch_*.txt

.PHONY: all race depends clean
