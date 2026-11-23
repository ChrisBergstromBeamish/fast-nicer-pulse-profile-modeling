# Setup

This page takes a fresh copy of the repository to the point where every check runs. Nothing here downloads anything for
you: the two public data sets the checks read are named with the addresses their own sources print, and you put them where
you like.

## 1. What is needed

- **The Python solution:** Python 3.12 and NumPy. Nothing else. (Written and checked with Python 3.12.10 and NumPy 2.5.4.)
- **The compiled solution:** a C compiler (C99) with the standard C and math libraries, and GNU make. Built and checked with
  gcc 13.3 on Ubuntu 24.04. On Windows the compiled solution is built and run on the Linux side, in a distribution under WSL
  (Windows Subsystem for Linux); the Python solution runs on the Windows side. A Python with NumPy is needed once at build
  time, to write the compiled program's table file (section 5 says how the two sides find each other).

## 2. The public reference material

The checks read two public data sets. Neither is copied into this repository (their licences and the choice of ours are
noted in `LICENCE_NOTES.md`).

1. **The data files of X-PSI's quick example.** X-PSI is the public code this work is compared with. Its source is at
   `https://github.com/xpsi-group/xpsi` (the address printed in its `pyproject.toml` and `README.rst`). The checks were
   written against the source tree at commit `bb5f461e605b22c6251b155b27a914a558757f93`, and only read its folder
   `examples/examples_fast/Data` (four text files, of which three are read: the example's counts, its background, its
   expected counts). Get that
   tree, for example with `git clone https://github.com/xpsi-group/xpsi` and `git checkout bb5f461e605b22c6251b155b27a914a558757f93`,
   and put it anywhere. X-PSI is not installed, not built and not run by anything here.
2. **The agreed reference curves and tables of Bogdanov et al. 2019**, "Constraining the Neutron Star Mass-Radius Relation
   and Dense Matter Equation of State with NICER. II. Emission from Hot Spots on a Rapidly Rotating Neutron Star",
   Astrophysical Journal Letters 887, L26 (article identifier ab5968; DOI 10.3847/2041-8213/ab5968). The journal distributes
   the article's supplementary data as the file set `apjlab5968`: 24 text files, described in its `ReadMe`. Put the folder
   that holds those 24 files anywhere.

## 3. Tell the programs where the reference material is

At the repository's top, copy `settings.cfg.example` to `settings.cfg` and fill in the two folders:

    xpsi = /where/you/put/xpsi            # the X-PSI source tree (its examples/examples_fast/Data is read)
    apjlab5968 = /where/you/put/apjlab5968   # the 24 files of the supplementary data

That is the only place the folders are given; no command below takes them as arguments. `settings.cfg` is this machine's and
is not committed (`.gitignore` lists it). On Windows, add one more line naming the WSL distribution in which the compiled
solution is built and run, for example `wsl = Ubuntu-24.04` (`wsl.exe --list` prints the names; without the line WSL's
default distribution is used). Folders may be written as Windows paths or as `/mnt/c/...` paths: the programs translate.

## 4. Check the Python solution

From the repository's top:

    python -I -B -X utf8 2-python-solution/check.py

It runs the ten quick checks and ends with the line `10 quick checks: 10 hold, 0 fail`, in about ten seconds on one
processor core of a laptop. Each line names what is compared, the difference found and the limit. The exit code is 0 when
every check holds. With `--slow` it also runs the slow checks (about ten minutes; `CLAIMS.md` lists them with their
times), and the last line then adds `; 20 slow checks: 20 hold, 0 fail`. Two slower comparisons are not in it and have their
own commands (`CLAIMS.md`). (The flags `-I -B -X utf8` make Python ignore any installed
packages and environment, write no cache files, and use UTF-8; every Python command in this repository is run that way.)

## 5. Build and check the compiled solution

On the Linux side, in the folder `3-compiled-solution`:

    make

This writes the table file `tables.bin` (the Python solution writes it: everything the compiled program needs that never
changes between workings-out), then builds the program in two forms, `nicer_fast` (the compiler's fast arithmetic, the form
timed so far) and `nicer_fast_strict` (strict arithmetic). To write the table file, `make` uses the first of `python3`,
`python` and `python.exe` that can import NumPy; on a Windows machine whose Linux side has no NumPy, `python.exe` is the
Windows Python, which WSL can run. If none is found, name one: `make PYTHON=/path/to/python`. `make depends` prints what each
built program is and which libraries it loads (`3-compiled-solution/README.md` records the answer).

Then, from the repository's top, on the side that has Python with NumPy:

    python -I -B -X utf8 3-compiled-solution/check.py

It compares both forms with the Python solution (the quick example, six stars at three settings, sixty random stars, 247 hard
made-up channels) and ends with `14 checks: 14 hold, 0 fail`, in under a minute. The exit code is 0 when every check holds,
and 2 if a program is not built.

**Windows notes.** Open the Linux side with `wsl.exe -d <distribution>`; a Windows folder `C:\a\b` is `/mnt/c/a/b` there.
The Python check commands run on the Windows side and call the compiled program through `wsl.exe` themselves. If `make`
warns about clock skew on a `/mnt/c` folder, the build is still complete; run `make` again to see it do nothing more.

## 6. The pages in this folder

- `DERIVATION.md`: what the programs compute, derived step by step from four named starting points, with the function that
  computes each formula and the check that tests it.
- `CLAIMS.md`: every number stated in this repository, in one table: the command that reproduces it, the value to expect,
  the allowed difference, the time it takes.
- `COVERAGE.md`: what the programs cover, where they refuse and how, what is approximate, what is not built, and the errors
  found so far.
- `PUBLIC_CODES.md`: the two public codes this work is set beside, X-PSI and the graphics-card code of Zhou and Huang, and
  why neither comparison is like for like yet.
- `LICENCE_NOTES.md`: the licences of the public sources read, and which of our files were written after reading which of
  theirs.
- `ROADMAP.md`: the sprints, one at a time.
