# Setup

This page takes a fresh copy of the repository to the point where every check runs. Nothing the checks need downloads
anything for you: the two public data sets the checks read are named with the addresses their own sources print, and you
put them where you like. (Only the install script of the comparison with the public code downloads, section 6.)

## 1. What is needed

- **The Python solution:** Python 3.12 and NumPy. Nothing else. (Written and checked with Python 3.12.10 and NumPy 2.5.4.)
- **The compiled solution:** a C compiler (C99) with the standard C and math libraries, and GNU make; the report of what was
  built (`make depends`) uses the common Linux tools `file` and `ldd`. Built and checked with gcc 13.3 on Ubuntu 24.04. On Windows the compiled solution is built and run on the Linux side, in a distribution under WSL
  (Windows Subsystem for Linux); the Python solution runs on the Windows side. A Python with NumPy is needed once at build
  time, to write the compiled program's table file (section 5 says how the two sides find each other).
  The batch mode (many stars at once) uses OpenMP, which gcc provides (`-fopenmp`, its library `libgomp`); its checks also
  need, on the Linux side, a `python3` with only its standard library, and gcc's race detector library (`libtsan`) with
  the tool `setarch`. All of these came with Ubuntu 24.04 and its gcc 13.3 package on the machine this was written on.

## 2. The public reference material

The checks read two public data sets. Neither is copied into this repository (their licences and the choice of ours are
noted in `LICENCE_NOTES.md`).

1. **The data files of X-PSI's quick example.** X-PSI is the public code this work is compared with. Its source is at
   `https://github.com/xpsi-group/xpsi` (the address its `pyproject.toml` prints as its home page). The checks were
   written against the source tree at commit `bb5f461e605b22c6251b155b27a914a558757f93`, and only read its folder
   `examples/examples_fast/Data` (four text files, of which three are read: the example's counts, its background, its
   expected counts). Get that
   tree, for example with `git clone https://github.com/xpsi-group/xpsi` and `git checkout bb5f461e605b22c6251b155b27a914a558757f93`,
   and put it anywhere. X-PSI is NOT needed installed for the checks: none of them installs, builds or runs it. Only the
   comparison of times and bins with it (section 6) needs it installed.
2. **The agreed reference curves and tables of Bogdanov et al. 2019**, "Constraining the Neutron Star Mass-Radius Relation
   and Dense Matter Equation of State with NICER. II. Emission from Hot Spots on a Rapidly Rotating Neutron Star",
   Astrophysical Journal Letters 887, L26 (article identifier ab5968; DOI 10.3847/2041-8213/ab5968). The journal offers the
   article's supplementary data on the article's own page, which the DOI leads to, as the file set `apjlab5968`: an archive
   (its `ReadMe` calls it "the tar.gz file") of 24 text files and the `ReadMe`. That page is the only address the sources
   give; the `ReadMe` prints none. Unpack the archive; the 24 files (`SD1a_test_IM.txt` to `OS1l_test_GSFCS.txt` and the two
   `_CU` tables) must lie directly in one folder, and that folder may be anywhere.

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

It runs the ten quick checks and ends with a line that begins `10 quick checks: 10 hold, 0 fail` (the seconds it took
follow), in about ten seconds on one processor core of a laptop. Each line names what is compared, the difference found and the limit. The exit code is 0 when
every check holds. With `--slow` it also runs the slow checks (about ten minutes; `CLAIMS.md` lists them with their
times), and the last line then adds `; 26 slow checks: 26 hold, 0 fail`. Two slower comparisons are not in it and have their
own commands (`CLAIMS.md`). (The flags `-I -B -X utf8` make Python ignore any installed
packages and environment, write no cache files, and use UTF-8; every Python command in this repository is run that way.)

## 5. Build and check the compiled solution

On the Linux side, in the folder `3-compiled-solution`:

    make

This writes the table file `tables.bin` (the Python solution writes it: everything the compiled program needs that never
changes between workings-out), then builds the program in two forms, `nicer_fast` (the compiler's fast arithmetic, the form
timed so far) and `nicer_fast_strict` (strict arithmetic), both able to work out many stars at once on several threads
(the batch mode); `nicer_fast_no_openmp`, the fast form built without OpenMP, which does the same on one thread; and
`libnicer_fast.so`, the fast form as a library for Python (`3-compiled-solution/README.md`). To write the table file, `make` uses the first of `python3`,
`python` and `python.exe` that can import NumPy; on a Windows machine whose Linux side has no NumPy, `python.exe` is the
Windows Python, which WSL can run. If none is found, name one: `make PYTHON=/path/to/python`. `make depends` prints what each
built program is and which libraries it loads (`3-compiled-solution/README.md` states the answer).

Then, from the repository's top, on the side that has Python with NumPy:

    python -I -B -X utf8 3-compiled-solution/check.py

It compares both forms with the Python solution (the quick example, six stars at three settings, sixty random stars, 247 hard
made-up channels: sixteen checks), then checks the batch mode against the program itself: each form, and the form without
OpenMP, gives every star in a batch exactly the bits it gets alone, on 1 to 8 threads; the library gives the batch mode's
lines; and `make race` (which it runs) finds no race. It ends with a line that begins `24 checks: 24 hold, 0 fail`, in
about a minute. The exit code is 0 when every check holds, and 2 if a program is not built.

**Windows notes.** Open the Linux side with `wsl.exe -d <distribution>`; a Windows folder `C:\a\b` is `/mnt/c/a/b` there.
The Python check commands run on the Windows side and call the compiled program through `wsl.exe` themselves. If `make`
warns about clock skew on a `/mnt/c` folder, the build is still complete; run `make` again to see it do nothing more.

## 6. The comparison with the public code (needs X-PSI installed; the checks do not)

The checks above need nothing of X-PSI but its example's data files (section 2). The comparison of its values, bins and
times with ours, whose printouts of 23 November 2026 are in `measured-2026-11-23/`, needs X-PSI installed and run. Three
scripts in `2-python-solution/tools/` do that on Linux (on Windows, on the Linux side). Each works in the one folder named
as `BASE` at its top, the same in all three; change it in all three to put X-PSI in another folder that you can write to.

1. `bash 2-python-solution/tools/install_reference_code.sh` installs X-PSI into that folder, with its own Python
   (Miniforge) inside it, in the way X-PSI's `docs/source/install.rst` gives, at the commit named in section 2. Unlike
   everything else in this repository it downloads: Miniforge, X-PSI's source and the packages of X-PSI's environment. It
   ends with the line `INSTALL_DONE`.
2. The next two scripts write into the folder `2-python-solution/results`, which they make (`.gitignore` keeps their
   `.jsonl` and `.npz` files out of the repository). `bash 2-python-solution/tools/run_public_code_timing.sh` runs X-PSI's own tests
   of its two shipped likelihood values and of its likelihood step, then times X-PSI at its shipped settings and three
   finer ones, for both shipped checks, and saves its bins (`results/public_code_times.jsonl`).
   `bash 2-python-solution/tools/refine_public_code.sh` refines its settings one at a time (`results/public_code_refine.jsonl`).
3. On the side that has our Python with NumPy, from the repository's top (the compiled rows need `make`, section 5):

       python -I -B -X utf8 2-python-solution/tools/compare_with_public_code.py 2-python-solution/results/public_code_times.jsonl "on mains"

   prints X-PSI's values, bins and times beside ours, and the ratios of the times for the quick example. The text in
   quotes is printed as the first line: say in it whether the laptop is on battery or on mains power, which changes every
   time measured. `python -I -B -X utf8 2-python-solution/tools/ours_convergence_both_stars.py` needs no X-PSI: it prints
   ours at settings finer than ours converged, at both checks' stars.

## 7. The pages in this folder

- `DERIVATION.md`: what the programs compute, derived step by step from four named starting points, with the function that
  computes each formula and the check that tests it.
- `CLAIMS.md`: every number stated in this repository, in one table: the command that reproduces it, the value to expect,
  the allowed difference, the time it takes.
- `COVERAGE.md`: what the programs cover, where they refuse and how, what is approximate, what is not built, and the errors
  found so far.
- `PUBLIC_CODES.md`: the two public codes this work is set beside, X-PSI and the graphics-card code of Zhou and Huang:
  what was measured of X-PSI on the same laptop, and why the comparison with the graphics-card code is not like for like.
- `measured-2026-11-23/`: what the comparison with X-PSI printed on 23 November 2026, kept as printed (its `README.md`
  says what each file is).
- `LICENCE_NOTES.md`: the licences of the public sources read, and which of our files were written after reading which of
  theirs.
- `ROADMAP.md`: the sprints, one at a time.
