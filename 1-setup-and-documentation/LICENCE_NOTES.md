# Licence notes: the licences of the public sources read, and which of our files were written after reading which of theirs

Facts only. This page says what the licences of the public sources are, as their own files state them, and which files of
this repository were written after reading which files of theirs. It gives no advice. The author (Chris Bergstrom) has chosen the
GNU General Public License, version 3, for this repository: the file `LICENSE` at its top.

## 1. X-PSI

The public code this work is compared with. The commit read: `bb5f461e605b22c6251b155b27a914a558757f93` (its
`pyproject.toml`: `name = "xpsi"`, `version = "3.5.0b1"`).

What its own files state:

- `pyproject.toml`, line 21: `license = "GPL-3.0-or-later"`; line 22: `license-files = ["LICENSE", "COPYING"]`.
- `LICENSE`: "Copyright (c) 2022-2026 the X-PSI core team", then "This program is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either
  version 3 of the License, or (at your option) any later version." It adds an "Additional permission under GNU GPL version
  3 section 7" for a work that links or combines the program with MultiNest.
- `COPYING`: the text of the "GNU GENERAL PUBLIC LICENSE", "Version 3, 29 June 2007".

So X-PSI's licence, in its own words, is the GNU General Public License, version 3 or (at your option) any later version.

What of it is in this repository: no file of X-PSI, and one stretch of lines. `2-python-solution/tools/time_public_code.py`,
the stopwatch of the comparison with the public code, builds the example's objects with lines copied from the example's
`examples/examples_fast/Modules/main.py`, between two comment lines that say so, so that it times what the example computes.
X-PSI is not imported, built or run by our programs or by the checks. It is installed, built and run only by the scripts
of the comparison (`install_reference_code.sh`, `run_public_code_timing.sh`, `refine_public_code.sh`, and
`time_public_code.py`, which imports it inside X-PSI's own environment), in a folder of its own, and only when a user runs
them (`README.md` of this folder, section 6). The checks read the text files of its shipped example's data folder (`examples/examples_fast/Data`) as numbers at run
time, from wherever the user has put the source tree (`README.md` of this folder, sections 2 and 3). Which of our files were
written after reading which of its source files is in section 5.

## 2. The GNU Scientific Library (GSL)

Our likelihood step follows the algorithm of GSL's Akima-periodic interpolation. `public_likelihood.py` says of its
`akima_periodic(x, y)`: "the interpolant of the GSL library called "akima-periodic" (algorithm read from GSL's akima.c)".
The algorithm was written again in Python (`AkimaPeriodic` in `public_likelihood.py`); no file of GSL is copied here, and
GSL is not linked or called by anything here. It is used in the checks of our likelihood step against the public code's
shipped unit values and in the "emulated" path of `quick_example.py`, which applies the public code's own interpolation
steps to our curve; the timed working-out does not use it.

GSL is not in the reference material on this machine: X-PSI's source tree only declares the GSL routines it calls, in
`xpsi/include/GSL.pxd`, and the file `akima.c` is not in that tree. So GSL's licence file was not re-read for this page.
GSL is distributed under the GNU General Public License, version 3 or later, as its source files state. That is stated here
from knowledge, and is still to be confirmed against GSL's own `COPYING` file.

## 3. The supplementary data of Bogdanov et al. 2019

"Constraining the Neutron Star Mass-Radius Relation and Dense Matter Equation of State With NICER. II. Emission From Hot
Spots on a Rapidly Rotating Neutron Star", Astrophysical Journal Letters 887, L26; the journal's file set `apjlab5968`,
24 text files with a `ReadMe`. The checks read these files as numbers at run time (the agreed reference curves, the table of
the deflection angle and lensing factor, the table of the travel time delay); none is copied here.

Its `ReadMe` states no licence. It describes the files and says: "These files are being provided to the community to
facilitate testing of other independently developed ray-tracing codes." That is the ReadMe's only statement about their
use.

## 4. The graphics-card code of Zhou and Huang (`gpu_ppm`)

The folder `gpu_ppm_3aaa0a8` of the reference material. Its `README.md` names it "gpu_ppm: GPU-Accelerated Pulsar / Neutron
Star Pulse Profile Modeling" and refers to the paper arXiv:2510.07764. Its `LICENSE` file is the MIT License:
"Copyright (c) 2025 Tianzhe Zhou and Chun Huang", followed by the MIT permission notice. The README also says that the
archive `ext.tar.xz` in the folder redistributes third-party data (atmosphere tables, the supplementary data of Bogdanov
et al. 2019, NICER response files and others) "verbatim", each under its original licence as given on its source page; that
archive was not opened.

Of this code, only `README.md` and `LICENSE` were read. Nothing of it is used, read or run by anything in this repository.

## 5. Which of our files were written after reading which of their files

"Written after reading" means: the author read the named source file, or the specification's verbatim quotes of it, and
then wrote our file in our own words and code; no line of theirs is copied. The specification named below,
`XPSI_MODEL_SPEC.md`, is a document written from X-PSI's source by reading it only (its own words: "Nothing in the source
tree was run, built, installed or imported"); it is not part of the repository.

| our file | written after reading | what was taken as the target or the definition |
|---|---|---|
| `2-python-solution/code/public_likelihood.py` | X-PSI's likelihood sources, as the specification quotes them in its sections 8.3, 8.5 and 9: `xpsi/Signal.py`, `xpsi/tools/energy_integrator.pyx`, `xpsi/likelihoods/default_background_marginalisation.pyx`, `examples/examples_fast/Modules/CustomSignal.py`; and GSL's `akima.c` | the definition of the number to be reproduced: the interpolation in energy and in phase, the search for the background's best value and limits, the integral per channel, the early exit and the failure value |
| `2-python-solution/code/quick_example.py` | `examples/examples_fast/Modules/main.py`, `xpsi/global_imports.py`, `xpsi/surface_radiation_field/hot_BB.pyx` | the example's parameter vector, its data files' names, the way its instrument edges are made, its exposure and effective area, its constants, its expected likelihood and tolerance |
| `2-python-solution/code/fast_kernel.py`, `2-python-solution/code/plain_star.py` | the coefficients of the star's flattened shape, which the specification quotes from `xpsi/surface_radiation_field/effective_gravity_universal.pyx` and `xpsi/cellmesh/mesh_tools.pyx` | the two coefficients −0.788 and 1.030 of the shape formula, an input of the example flagged in `COVERAGE.md` section 1.4; everything else in these two files is ours |
| `2-python-solution/check.py`, `2-python-solution/code/test_public_likelihood.py` | `xpsi/tests/test_background_marginalisation.py` | the unit test's inputs (its counts, its signal at 64 phases, its exposure and phase shift) and its two expected values of the likelihood step |
| `2-python-solution/tools/second_check_point.py`, `2-python-solution/tools/second_check_point_emulated.py` | `xpsi/tests/test_likelihood.py` | the second check's star, its setting (no upper limit on the background, 128 energies), its expected value and tolerance |
| `2-python-solution/tools/quick_example_counts_check.py` | the example's shipped file of expected counts, `xpsi_good_expected_hreadable.dat` (read at run time, not copied) | the target of the bin-by-bin comparison of counts |
| `2-python-solution/tools/time_public_code.py` | `examples/examples_fast/Modules/main.py`; `xpsi/tests/test_likelihood.py` | lines copied from the example's `main.py` (its data, its instrument, the limits of its background, its star, its hot region and its likelihood object: the stretch between the comment lines "copied from main.py" and "end of the copied part"), and the second check's star and shipped value from the test file. The script runs only inside X-PSI's own environment, for the comparison |
| `3-compiled-solution/nicer_fast.c` | our own `public_likelihood.py` and `fast_kernel.py`, of which it is a port ("A line-by-line port of our own Python programs"); so, at one remove, the same sources as those two files | the same definitions, at one remove |
| `XPSI_MODEL_SPEC.md` (the specification; not in the repository) | the X-PSI files it quotes verbatim, listed below | the description of what one likelihood evaluation of the example computes |

The X-PSI files the specification quotes verbatim (its headers of the form `file : lines`), at commit
`bb5f461e605b22c6251b155b27a914a558757f93`:

- `examples/examples_fast/Modules/main.py`
- `examples/examples_fast/Modules/CustomInstrument.py`
- `examples/examples_fast/Modules/CustomPhotosphere.py`
- `examples/examples_fast/Modules/CustomSignal.py`
- `xpsi/__init__.py`
- `xpsi/HotRegion.py`
- `xpsi/Likelihood.py`
- `xpsi/ParameterSubspace.py`
- `xpsi/Photosphere.py`
- `xpsi/Signal.py`
- `xpsi/Spacetime.py`
- `xpsi/Star.py`
- `xpsi/global_imports.py`
- `xpsi/cellmesh/integrator_for_azimuthal_invariance.pyx`
- `xpsi/cellmesh/mesh.pyx`
- `xpsi/cellmesh/mesh_tools.pyx`
- `xpsi/cellmesh/rays.pyx`
- `xpsi/likelihoods/default_background_marginalisation.pyx`
- `xpsi/surface_radiation_field/effective_gravity_universal.pyx`
- `xpsi/surface_radiation_field/hot_BB.pyx`
- `xpsi/surface_radiation_field/hot_wrapper.pyx`
- `xpsi/tools/compute_expected_counts.pyx`
- `xpsi/tools/core.pyx`
- `xpsi/tools/energy_integrator.pyx`

The specification also names, without quoting them, `xpsi/include/GSL.pxd` (the declarations of the GSL routines X-PSI
calls) and the plain-text file `.git/packed-refs` (the source of the commit identifier).

## 6. Statements

- No file of X-PSI, of GSL or of the article's supplementary data is copied into this repository. The two data sets the
  checks read are read at run time from the folders the user names in `settings.cfg`.
- X-PSI is not imported, built or run by our programs or by the checks. The scripts of the comparison with it install,
  build and run it in a folder of its own when a user runs them; one of them, `time_public_code.py`, imports it there and
  holds lines copied from its example's `main.py` (sections 1 and 5).
- Nothing of the graphics-card code of Zhou and Huang is used.
- The licence of this repository: the GNU General Public License, version 3, chosen by the author (Chris Bergstrom) on
  23 November 2026; the file `LICENSE` at the repository's top holds its text.
