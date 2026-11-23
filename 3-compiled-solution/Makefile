# The compiled solution: builds the table file and the program in two forms. Run on the Linux side, in this folder:
#     make            writes tables.bin (once), builds nicer_fast (the compiler's fast arithmetic, the form timed so far) and
#                     nicer_fast_strict (strict arithmetic: no -ffast-math, and no fused multiply-add either)
#     make depends    says what each built program is and which libraries it loads (file, ldd)
#     make clean      removes what make made
# Then, from the repository's top, on the side that has Python with NumPy:  python -I -B -X utf8 3-compiled-solution/check.py
# tables.bin is written by the Python solution (../2-python-solution/tools/export_for_c.py) and needs a Python with NumPy. The
# first of python3, python and python.exe that can import NumPy is used; on the Linux side of a Windows machine whose Linux
# Python has no NumPy, python.exe is the Windows Python, reached through WSL. Or name one:  make PYTHON=/path/to/python
# The public example's data folder comes from ../settings.cfg (1-setup-and-documentation/README.md).
CC = gcc
CFLAGS_FAST = -O3 -march=native -ffast-math -Wall -Wextra
CFLAGS_STRICT = -O3 -march=native -ffp-contract=off -Wall -Wextra
PYTHON ?= $(shell for p in python3 python python.exe; do if $$p -I -c "import numpy" >/dev/null 2>&1; then echo $$p; break; fi; done)

all: tables.bin nicer_fast nicer_fast_strict

nicer_fast: nicer_fast.c
	$(CC) $(CFLAGS_FAST) -o $@ $< -lm

nicer_fast_strict: nicer_fast.c
	$(CC) $(CFLAGS_STRICT) -o $@ $< -lm

tables.bin: ../2-python-solution/tools/export_for_c.py ../2-python-solution/code/fast_kernel.py ../2-python-solution/code/quick_example.py ../2-python-solution/code/public_likelihood.py ../2-python-solution/code/consts.py
	@test -n "$(PYTHON)" || { echo "no Python with NumPy was found (tried python3, python, python.exe): run  make PYTHON=<a python with NumPy>" >&2; exit 1; }
	$(PYTHON) -I -B -X utf8 ../2-python-solution/tools/export_for_c.py $@

depends: nicer_fast nicer_fast_strict
	file nicer_fast nicer_fast_strict
	ldd nicer_fast
	ldd nicer_fast_strict

clean:
	rm -f nicer_fast nicer_fast_strict tables.bin hard_channels*.bin hard_channels*.bin.out dump_*.bin

.PHONY: all depends clean
