# Coverage: what the programs cover, where they refuse, what is approximate, what is not built, and the errors found so far

This page is written for a reader who wants to check rather than trust. Each point is put in plain words first and then as
the exact statement in the program, or the exact measurement, it rests on. The programs are three: the fast program in
Python (`2-python-solution/code/fast_kernel.py`, function `spot_counts`, with the public likelihood's definition in
`public_likelihood.py` and the public example's inputs in `quick_example.py`); the plain way (`plain_star.py` for a
flattened star, `plain_sphere.py` for a round one), which is slow and exists to check the fast program bin by bin; and the
compiled copy (`3-compiled-solution/nicer_fast.c`), the fast program written out again in C. Every measured number on this
page has a row in `CLAIMS.md` with the command that reproduces it; where a script is named in brackets it is in
`2-python-solution/tools/`, and `2-python-solution/check.py` and `3-compiled-solution/check.py` run the checks in one command
each.

Agreement with the public code (X-PSI) is not a result about the physics and never a proof of it. `DERIVATION.md` says from
which starting points the programs compute light's path, and why this project cannot tell those starting points apart from
the public code's. This page compares; the only things of the public code's in our programs are the inputs of its example,
which are flagged where they appear.

## 1. What the programs cover

**In plain words.** One circular hot spot of one temperature on a spinning star, seen from any direction, with one path of
light from every point of the spot to the observer; a blackbody spectrum; and the instrument, data and constants of the
public code's shipped quick example. The star may be round or flattened by its spin; the flattened shape is the public
example's own formula, used as an input and flagged below.

### 1.1 The star's compactness, L/R

L is the length 2GM/c² for the star's mass and R its radius. The programs form the ratio at the equator and at the pole
(`u_eq`, `u_pole` in `spot_counts`), and both must lie inside the range.

| program | range of L/R | the statement in the code |
|---|---|---|
| Python, with the path integrals themselves (`bend`) | 0.12 to 0.62 | `fast_kernel.py`: `U_RANGE = (0.12, 0.62)    # L/R covered by the program`; refused outside it (section 2) |
| Python, with the approximation of the path integrals (`_path_fit`; the default, `fit=True`) | 0.15 to 0.55, and every ring's tilt with the sine of its angle below 0.08 | `fast_kernel.py`: `FIT_U = (0.15, 0.55)`, `FIT_TILT = 0.08`; when either condition fails the Python uses the integrals instead, so the star is still worked out |
| compiled | the approximation's range only | `nicer_fast.c`, `spot_counts`: `FIT_U0 < u_min && u_max < FIT_U1 && fabs(o2) < 0.5` and `h_max < FIT_TILT`, the values read from the table file that `tools/export_for_c.py` writes from the Python's `FIT_U` and `FIT_TILT`; refused with a code otherwise (section 2) |

The approximation of the path integrals is one two-variable Chebyshev series for the swept angle and the delay, made once
from the programs' own integrals (`bend`). What it costs in accuracy is in section 3.

### 1.2 One path of light

Every point of the spot sends light to the observer along one path. Where a second path around the star is possible, the
programs stop instead of giving a number (section 2 has the test and its message). The `TwoPaths` docstring in
`fast_kernel.py`: "The star is so compact (or so flattened) that light from behind reaches the observer by a second path
around it. The program adds one path only and refuses such a star instead of returning a wrong number."

### 1.3 The spot, the observer, the spectrum

- One circular spot of angular radius `rho` centred at colatitude `theta_s` (`spot_rings`: "a circular spot"); it may cover
  a pole (the latitude range is cut where whole rings begin). One temperature `kT` for the whole spot.
- The observer at any inclination from 0 to π, including exactly on the spin axis (`aligned`, `B < ALIGN_B` with
  `ALIGN_B = 1.0e-6`: "a ring is taken as seen from along the axis").
- A blackbody spectrum: `bb_photons(E, kT) = BB_CONST * E * E / expm1(E / kT)`, the same in every direction (the `beaming`
  argument is `None` in every use here).
- Geometries checked against the plain way and against the agreed reference curves of Bogdanov et al. 2019 (the measured
  results, `kernel_against_reference_curves.py`, `kernel_against_plain_way.py`, `kernel_against_plain_star_many.py`,
  `kernel_against_plain_star.py`):
  - the sixteen agreed reference curves, directly: twelve extended spots within 1.2e-5 to 2.8e-4 of the curve's median
    (the article's target is 1e-3); four point-like spots within 2.4e-3 to 7.7e-3 at the edges of eclipses, which is above
    that target there (the project's record states that the article reports the same of its own codes; the article itself
    is not in the reference material on this machine, so that statement was not re-read for this page);
  - round stars, 36 geometries (observers from 0° to 180°, also exactly on the axis; spots in both halves and over the
    poles): 34 within 2e-4 of the largest bin (most near 1e-5); two at 2.4e-4 and 3.2e-4, where only a sliver of the spot
    is in view, and there the difference shrinks as the plain way's patches are made smaller (to 1.1e-4 and 2.1e-4 at
    four times as many patches: the project's record; the committed script runs one patch setting);
  - flattened stars, 24 geometries (200 to 716 Hz, observers at any angle and near the axis; the committed script
    `kernel_against_plain_star_many.py`, about fifteen minutes, run on its own): 22 worked out, all within 1.9e-5 of the
    largest bin; 1 refused by the fast program itself (a second path of light is possible there); 1 could not be checked,
    because the plain way refused it (its table of directions reaches further inward than any visible ray). The project's
    record speaks of 29 geometries with 27 within 3.2e-5; the committed script holds 24, and its numbers are the ones
    stated here;
  - 60 stars drawn at random from the whole range the public example's prior allows: both the Python and the compiled
    program ran without stopping on every one (`prior_draws_c_and_python.py`). "The public example's prior stays below
    0.4725 and is never refused."

### 1.4 The public example's instrument, data and constants, used as they are

The checks read the data files of X-PSI's quick example (`examples/examples_fast/Data`) as numbers; nothing of X-PSI is
imported or run (`quick_example.py`: "X-PSI is not imported and not run. Its data files are read as numbers; its constants
are used as the example's inputs").

- The instrument: 291 input intervals from 0.15 to 1.605 keV, made the way the example's `main.py` makes them
  (`example_instrument_edges`: repeated addition of 0.005 from 0.1 and 0.105, rows 10 to 300 kept); one constant effective
  area, `AREA = 1800.0`, which our program folds with the energy integrals into one fixed matrix.
- The data: 291 channels by 32 phase bins, an exposure of 1000 s, 578,528 counts; the background in each channel limited
  to 0.9 to 1.1 times the shipped background (`load_data`: `support = [0.9 * bg / T_EXP, 1.1 * bg / T_EXP]`).
- The star, in the example's own parameter names: `mass` 1.4 (solar masses), `radius` 12 (km), `distance` 1.0 (kpc),
  `cos_inclination` cos 60°, `p__phase_shift` 0.0, `p__super_colatitude` 70° in radians, `p__super_radius` 0.75 rad, and
  the eighth entry of its vector `p`, 6.7, the log10 of the spot's temperature in kelvin; the spin 314 Hz.
- The public code's constants, used as the example's inputs (`quick_example.py`, "inputs of the example, not ours"):
  `X_GM_M = 0.5 * 2.95325024e3` (GM_sun/c² in metres), `X_KPC_M = 3.08567758e19`,
  `X_KB_OVER_KEV = 1.38064852e-23 / 1.60217662e-16`, `X_BB = 5.040366110812353e22 * 1.0e-7 / 1.60217662e-16`. The star's
  mass is turned into the length L with the public code's value of GM_sun/c² (`star_args`), so that the comparison is of
  the same star.
- **The star's flattened shape: an input of the public example, flagged.** The programs take
  `R(θ) = R_eq (1 + o2 cos²θ)` with `o2 = (Ω² R_eq³ / GM) (−0.788 + 1.030 GM / (R_eq c²))` (`fast_kernel.py`,
  `plain_star.py`, `nicer_fast.c`, the same line in each). This formula and its two coefficients are the public example's
  definition of its star, written to match general relativity's star models; they are used only because the example
  defines its star so, and the derivation of light's path, the delay and the energy shifts does not use them (see
  `DERIVATION.md`). Nothing of ours rests on them.

### 1.5 In parallel: many stars at once

**In plain words.** The compiled program can work out many stars at the same time, one star on each processor core, and
each star comes out exactly as it does alone. One star is still worked out on one core.

"In parallel" here means only this: the batch mode of `nicer_fast` and `nicer_fast_strict`, and the library
`libnicer_fast.so`, work out many stars at once, one star per thread at any moment (threads are lines of work of one
program, each on a core; OpenMP, a feature of the compiler, hands the stars to them). For every star, ln L, the code and
every bin count come out the same, bit for bit, as when that star is worked out alone, on any number of threads and in any
order of the stars (`3-compiled-solution/check.py`; `3-compiled-solution/README.md`). It does not cover: parallel work inside
one star (the rings of one spot are worked out one after another); several stars inside one thread with the processor's
vector instructions; the Python solution, which works on one star at a time on one core; a graphics card; a build for
the Windows side (the compiled solution is built and run on the Linux side); and machines other than the one laptop on
which it was measured (`CLAIMS.md`, group I). The time on several cores is measured for our program only; X-PSI's time on
several threads is not measured, so no ratio against X-PSI is stated for more than one core.

## 2. Where the programs refuse, and how

**In plain words.** Where a star is outside what the programs cover, or where light could reach the observer by two paths,
or where a step of the arithmetic cannot be completed, the Python programs stop with a message and the compiled program
prints a code. Where the public code itself would leave early, our programs give −∞ (the Python's step-by-step form raises
an error instead, so that the case is seen). No wrong number is returned in any of these cases. Below, each message is
quoted as the fixed text in the source; `%.3f`, `%.4f`, `%d` and `%s` mark where the program prints its numbers.

### 2.1 `fast_kernel.py`, function `spot_counts`

| exception | message (fixed text) | when |
|---|---|---|
| `ValueError` | `L/R = %.3f to %.3f is outside the range the program covers (%.2f to %.2f)` | unless `U_RANGE[0] < min(u_eq, u_pole)` and `max(u_eq, u_pole) < U_RANGE[1]`, that is unless both ratios lie strictly inside 0.12 to 0.62; the two printed ratios are the lower and the higher of the two |
| `TwoPaths` | `a ring's tilt (sin %.3f) reaches the directions that do not come back out (L/R %.3f)` | for some ring, `x_floor < x_crit + 0.01`, where `x_floor = −(|sin tilt| + 0.02)` is the most inward direction in which a visible ray can leave ("no visible ray starts further inward than the tilt") and `x_crit = −sqrt(max(1 − 6.75 u² (1 − u), 0))` is the direction inward of which "a ray does not come back out" |
| `TwoPaths` | `light from behind the star reaches the observer by a second path: at latitude %.3f rad a visible ray can sweep %.3f rad, and the far point of that ring is %.3f rad from the observer's direction` | for some ring, `th_edge >= 2π − psi_far − 1e-9`, where `th_edge` is the swept angle of the ray that leaves edge-on (at `x = −|sin tilt|`) and `psi_far = arccos(cos i cos θ − sin i sin θ)` is the angle of the ring's farthest point from the observer's direction |
| `TwoPaths` | `a swept angle above pi on the first path (%.4f rad)` | a swept angle at one of the ring's points is not below π + 1e-9; the comment beside it: "cannot happen after the test above; kept as a guard" |
| `ArithmeticError` | `an end direction of a ring could not be bracketed (latitudes %s rad)` | after the secant steps a needed end direction of a ring still has a residual of `END_TOL = 1.0e-12` or more, and the function has the same sign at both ends of that direction's bracket, so no bisection can find it |

`TwoPaths` is a subclass of `NotImplementedError`. The test behind the second `TwoPaths` message is the one the results
describe, and it is quoted here exactly: "The test: no visible ray sweeps more than the one that leaves edge-on; if that
sweep reaches 2π minus the angle of a ring's farthest point, a second path is possible and the program stops with a message
instead of giving a number. For a round star this can first happen at L/R = 0.568, and only for spots that reach the far
side. The public example's prior stays below 0.4725 and is never refused." The comment block "one path of light only" in
`spot_counts` gives the reasoning: "A ray that leaves the surface has x > -|sin(tilt)|, and the swept angle grows as x
falls, so no visible ray sweeps more than Theta(u, -|sin(tilt)|). The farthest point of a ring's turn is at angle psi_far; a
second path to some point of the ring needs a sweep of at least 2 pi - psi_far. If the largest sweep is below that, there
is one path only; if not, the program refuses the star." The example that showed the earlier test to be incomplete, a round
star with L/R = 0.6 seen from its axis with a small spot near the far pole, is now refused; the same star with the spot on
the near side is still worked out (section 6, item 9).

### 2.2 `plain_star.py`, function `light_curve`

| exception | message (fixed text) | when |
|---|---|---|
| `fast_kernel.TwoPaths` | `the plain way also adds one path only` | any swept angle in the patch's table of directions reaches π (`not np.all(th < math.pi)`) |

### 2.3 `public_likelihood.py`

The likelihood is the public code's own definition, written out again (its docstring: "This file is the definition of the
number to be reproduced; it is not physics of ours and it holds no physics"). It has two forms that give the same number:
the step-by-step form (`marginal_loglike`, `marginal_loglike_channel`) and the fast form (`marginal_loglike_fast`, all
channels at once).

| where | exception or value | message (fixed text) or meaning | when |
|---|---|---|---|
| `marginal_loglike` (step by step) | `ValueError` | `channel %d: the public code would leave by its early exit (model far above the data)` | in some channel, `src − 20 sqrt(src) − (counts in the channel) > 0`, with `src` the model's counts in the channel (its rates summed over the bins, times the exposure over the number of bins) |
| `marginal_loglike_fast` | returns `−inf` | the public code's early exit; the comment: "the public code leaves here with a number near -1e90" | the same condition in any channel |
| `marginal_loglike_channel` | returns `−inf` | "the public code's failure exit (a number near -1e90)" | the logarithm of the integrand at the public code's best background is not finite: some bin's model rate plus background is zero or below, unless both that value and the bin's counts are within 1e-12 of zero |
| `marginal_loglike_channel` | `ArithmeticError` | `the reference integral did not settle` | the Gauss-Legendre panels of the channel's integral have been doubled to 16384 and two successive answers still do not agree to `rtol = 1.0e-12` |
| `marginal_loglike_channel` | `ArithmeticError` | `the reference integral came out zero or not finite` | the settled integral is not a positive finite number |
| `marginal_loglike_fast` | `ArithmeticError` | `a channel's integral came out zero or not finite in the fast form` | a channel's integral in the fast form is not finite |

The early exit and the failure exit are the public code's own rules, written again; they are not refusals of ours, and they
return what the public code returns at that point (a number near −1e90, written here as −∞). The step-by-step form raises
`ValueError` on the early exit instead, so that the case is seen when the two forms are compared.

### 2.4 `nicer_fast.c`

The function `spot_counts` "Returns 0, or a negative code: -1 outside the covered range, -2 a second path of light, -3 an
end direction not found." The main program then prints `refused: code %d` and exits with status 1.

| code | printed | when |
|---|---|---|
| −1 | `refused: code -1` | an input is not in its range: unless `mass > 0`, `Req_km > 0`, `D_cm > 0`, `kT > 0`, `0 ≤ incl ≤ π`, `0 < rho < π`, `0 < theta_s < π`; or unless `FIT_U0 < u_min && u_max < FIT_U1 && fabs(o2) < 0.5` (L/R inside 0.15 to 0.55 at both the equator and the pole, and the flattening `o2` below 0.5 in size); or unless `h_max < FIT_TILT` (the bound on the surface's tilt below 0.08); or `8 * NPP > MAXR` (more than 64 rings per panel: the arrays hold `MAXR = 512` rings, and the comment beside the test says "at most 7 panels of NPP rings") |
| −2 | `refused: code -2` | a second path of light is possible: `x_floor < x_crit + 0.01` for a ring; or the swept angle of the edge-on ray `>= TWO_PI − acos(A − B) − 1e-9` (the same test as the Python's); or, as a guard, a swept angle at a point not below π + 1e-9 |
| −3 | `refused: code -3` | the search for an end direction of a ring (the direction of its nearest point, of its farthest point, or the edge-on direction between them) did not find it |

`HUGE_NEG` is defined as `(-1.0e300)` and "stands for \"minus infinity\" (the public code's failure value is near -1e90)".
The main program prints `ln L = -inf` when the likelihood is below −1e299. That happens:

- at the public code's early exit, "the model far above the data" in some channel (the same condition as the Python's);
- in a channel sent to the reference form, when the logarithm of the integrand is not finite, when its integral did not
  settle in 16384 panels, or when it came out zero (the comment: "did not settle, or zero: refused (the Python raises
  here)");
- in the fast form, when a channel's integral is zero ("an integral of zero: the public code's failure exit").

So the compiled program answers `ln L = -inf` in two cases where the Python raises `ArithmeticError`. The checks compare
the two programs only where both return a finite likelihood (`3-compiled-solution/check.py`: "sixty prior draws with a
finite likelihood in both").

## 3. What is approximate

**In plain words.** Two steps of the fast program are exact identities (counting the photons when they leave removes the
lensing factor; a whole stretch of the spot goes into a time bin exactly; `DERIVATION.md`). Everything around them is
arithmetic with a finite number of pieces, and each piece was varied on its own to measure what it costs. The measurements
are below, each with where it comes from.

### 3.1 The finite number of rings, points and energies

"Rings" are rings of latitude across the spot (`n_rings` per panel of the latitude range); "points" are directions in
which the light leaves, per ring and half-turn (`n`); "energies" are the energies at which the spectrum is worked out
before it is integrated over each channel (`nE`). How the quick example's likelihood moves with them (`quick_example.py`
prints these rows; `check.py` holds the converged value against the shipped one and the working setting against the
converged value):

| how | ln L | minus X-PSI's shipped value |
|---|---|---|
| X-PSI's shipped value (its own tolerance: ±0.3255) | −32552.616905 | 0 |
| ours, 48 rings per panel, 32 points, 48 energies | −32552.710363 | −0.0935 |
| ours, 24 rings per panel, 24 points, 32 energies | −32552.709070 | −0.0922 |
| ours, 12 rings per panel, 16 points, 24 energies | −32552.711306 | −0.0944 |
| ours, 10 rings per panel, 16 points, 24 energies (the working setting) | −32552.719298 | −0.1024 |
| ours, 8 rings per panel, 12 points, 16 energies | −32552.774590 | −0.1577 |
| ours, 6 rings per panel, 12 points, 16 energies | −32553.007080 | −0.3902 (outside the tolerance) |
| our second program, the plain way, 48 × 192 patches | −32552.7071 | −0.0902 |
| our second program, the plain way, 96 × 384 patches | −32552.7095 | −0.0926 |
| X-PSI's two interpolation steps applied to our curve | −32552.669 | −0.052 |

Our value is steady to about 0.0022 from 12 rings per panel upward. The two programs of ours agree to 0.001 in the score
and to 4 parts in a million of the mean bin in every bin. They add the light up in different ways and share the path
integrals and the star's formulas (shape, tilt, area, energy shift, Doppler factor, the factors of the count formula), so
their agreement checks the adding-up, not those shared formulas; the agreed reference curves (section 1.3) and the counts
file (section 5) check those.

### 3.2 The even grid of 129 points that locates arrival times

The point where a stretch's arrival time meets a bin edge is found by a table on an even grid in the Chebyshev variable
and is not refined (`ND = 129`). The docstring of `fast_kernel.py`: "for the exact integrand an error in that point changes
K only at second order (the summand vanishes there); for the interpolated integrand there is also a first-order term, the
interpolation's residual at that point times the error". Measured by varying the grid alone from 65 to 2049 points
(`arrival_grid_check.py`, part of `check.py --slow`): against the 2049-point grid, the working 129 points change the quick
example's ln L at the working setting by 1.8e-6 and the bins by 3.5e-9 of the mean bin; 65 points change ln L by 1.3e-5 and
the bins by 2.4e-8. (The docstring gives 2e-6 and 4e-9 at 129 points; the project's record gives 9e-6 at 65 points, where
the committed check gives 1.3e-5.)

### 3.3 The path integrals and their approximation

- The path integrals themselves (`bend`) are Gauss-Legendre sums with `n_gl` points (32 in `spot_counts`, 96 where the
  tables and the approximation are built). The docstring of `spot_counts`, measured on the quick example and on the second
  check star: "32 points change ln L by 7e-12 and 1e-6 against 96 points; 3 steps by 7e-12 and 2e-6 against 8 steps" (the
  steps are the secant steps for each ring's end directions, `n_sec = 3`, whose residuals are then checked to `END_TOL`).
  These two measurements stand in the docstring only; no committed command reproduces them.
- The approximation of the path integrals (`_path_fit`): one two-variable Chebyshev series (32 by 72 terms) in u = L/R
  and z = sqrt(1 − x) for the swept angle over z and the delay over z², built once from `bend` with 96 points, and used
  for every star inside `FIT_U` with every ring's tilt below `FIT_TILT`. Measured against `bend` with 96 points on 6000
  random directions (`path_approximation_check.py`, part of `check.py --slow`): largest difference 8.5e-13 rad in the
  swept angle and 1.6e-12 R/c in the delay (the comment in `fast_kernel.py` gives 9e-13 rad and 1.3e-12 R/c from another
  sample of directions). The integrals themselves with 96 points against 160 points, on the same directions: 8.4e-14 rad
  and 4.8e-14 R/c.

### 3.4 The latitude cuts

The latitude range of the spot is cut where a ring's visibility changes kind (up to four latitudes), because a ring's
contribution has a corner there (`visibility_breaks`). Where the cuts are placed is approximate:

- From tables where the surface's tilt is small: the swept angle near the edge-on direction comes from tables of u alone
  "to second order in the tilt (good to about 2e-5 rad for a tilt of 0.03, to about 0.02 rad for a tilt of 0.2)". The
  programs use the tables alone when the bound on the tilt `h_max` is 0.05 or less ("below a tilt of 0.05 the tables place
  the cuts to better than 1e-4 rad"). The measured placing, from the review that asked for it: "from tables for slow and
  moderate spins (within 1.3e-4 rad of the exact place at 314 Hz and below)".
- On the exact function where the tilt is larger (`h_max > 0.05`): the contacts are moved onto the swept angle itself by
  further secant steps, "within 1e-8 rad" of the exact place, "solved to convergence".
- What a missing or misplaced cut costs was measured when all four kinds were placed: spots that are almost wholly hidden
  went from 1% error to 1e-4.

The three measurements of this section (1.3e-4 rad, 1e-8 rad, 1% to 1e-4) are the project's record of the work that placed
the cuts; no committed command reproduces them on their own. What the cuts leave in the bins is covered by the comparisons
with the plain way (section 1.3), which have committed commands.

### 3.5 Rings whose visible part is in two pieces

The fast program takes the visible part of each ring's turn as one stretch. Quoted exactly from the results: "The visible
part of each ring's turn is taken as one stretch. Checked on 18,000 rings of 300 stars in the public example's prior and on
3,000 latitudes of its star: no exception. At 500 Hz and above, a few rings per thousand have a visible part in two pieces
(largest μ there 0.04); the program then sets the hidden middle to zero at its points, which is approximate, and it has no
bound on the error this leaves in a bin. In the geometry test against the plain way, the three such stars that were tried
came out within 1.4e-5."

The committed test (`one_visible_stretch_check.py`, part of `check.py --slow`) draws 400 stars from a wider range than the
public example's prior (1.0 to 2.2 solar masses, 9.5 to 15 km, 100 to 716 Hz) and finds, of 14,968 rings, 13 whose visible
stretch does not start at the ring's nearest point. The largest share of a ring's light in such a stretch is the whole of it,
for a ring at 164° on a star of 14.3 km spinning at 716 Hz, whose surface there has a tilt with sine 0.21: beyond the 0.08
the approximation of the path integrals accepts, so the compiled program refuses that star, and the Python program works it
out with the integrals and takes that ring as one stretch starting at its nearest point, where it is hidden; the largest μ
on that ring is 0.0175. The record's statement about the public example's prior alone ("no exception") has no committed
command that reproduces it as stated; the committed test's numbers are the ones a reader can repeat.

### 3.6 The likelihood step

The likelihood's integral over the background rate of each channel is the public code's own definition; our fast form
computes it with one Gauss-Legendre rule over the part of the range where the integrand is above e^−36 of its largest
value, with 10 to 64 points chosen from how far the logarithm falls, and the cut points checked against the true logarithm
and moved outward until the fall is at least 32 (`DROP = 36.0`, `CUT_FALL = 32.0`). Measured against the step-by-step form
on 247 hard made-up channels (few counts, source far above or below the data, wide and narrow limits): 6e-11 per channel
(`fast_likelihood_stress.py`); the source's own words: "a measurement on those channels, not a bound for every channel".
The step-by-step form doubles its panels until two answers agree to 1e-12. Against the public code's two shipped values for
this step alone (tolerance 1e-8): −4521.763801535126 against −4521.763801536406 (relative difference 2.8e-13) with the
background limited to [0, 1e-11], and −2772.718471175409 against −2772.718471175409 (0) with no upper limit.

## 4. What is not built

**In plain words.** The programs stop at the edge of what they cover instead of going on. The list below is what lies past
that edge, with the reason it was left there.

- **The second path of light around very compact stars.** From L/R = 0.568 for spots on the far side (for a round star;
  earlier for a flattened one), light from behind the star can reach the observer by a second path. The fast program
  refuses every such star (section 2); working the second path out is not built. The public example's prior never reaches
  it.
- **A tabulated hot-gas (atmosphere) spectrum in place of the blackbody.** The programs hold the blackbody only
  (`bb_photons`); the public code's tabulated atmospheres are not in them, and the comparison so far is on the shipped
  example, which is a blackbody.
- **A test of whether anything stands in the light's way.** Quoted exactly: "Whether a ray, once it has left the surface,
  is caught again by the bulging equator is not tested (the public code does not test it either)." The programs also take
  it as a named premise that every photon that leaves in an open direction arrives (nothing absorbs it on the way).
- **Rings whose visible part is in two pieces.** Section 3.5: at 500 Hz and above a few rings per thousand; the hidden
  middle is set to zero at the ring's points, without a bound on the error. Handling the two pieces as two stretches is not
  built.
- **The background step without an upper limit is slow in Python.** With no upper limit on the background (the setting of
  the public code's second check) each channel's integral spans a wide range: the Python likelihood step alone then took
  5 to 8 ms, and the compiled copy 0.9 to 1.1 ms for the whole working-out, on 22 November 2026 (`CLAIMS.md`, row F6, has
  the slower repeat of 23 November); its values there equal the Python's to 1e-8. Fixed points for that case are not built.
- **Several hot spots, or shapes other than a circle.** `spot_counts` works out one circular spot per call.
- **Parallel work inside one star, several stars per thread with vector instructions, a graphics card.** Section 1.5: the
  compiled program works out many stars at once, one per thread; nothing else is built in parallel. The graphics card is
  the roadmap's next item.
- **Listed in the fourth review as not acted on, with the reasons given there** ("listed here so that they are not
  forgotten"):
  - exact handling of a ray that leaves exactly radially or exactly along the surface in the integrals (the approximation
    of the path integrals returns zero at the radial end; the tangent values come from a ray 1e-12 away);
  - a certificate for the inward turning point inside the integrals (the kernel only asks inside its guarded range);
  - the public code's own equality test and its handling of several components, energies beyond its grid and failure
    values (none is active in the shipped example).
- **The public code beyond its two shipped checks.** X-PSI was installed and run on this laptop on 23 November 2026
  (`PUBLIC_CODES.md`, section 1.2), at the quick example's star and the second check's, at its shipped settings and 31
  finer ones (`CLAIMS.md`, group H). Not done: a comparison with it at other stars, with several spots or a tabulated
  atmosphere, or on another machine; and a time of it at the second check's star (section 6, item 13).

## 5. What the agreement with the public code does not show

**In plain words.** Our likelihood for the shipped quick example passes the public code's own check, and our two programs
agree with each other bin by bin. That is a result about a method of calculation. With the public code run on this laptop
(23 November 2026), its likelihood step gives ours to the last printed digit, so every difference in the score lies in the
counts; and as its settings are refined its counts move toward ours, per phase bin with its mesh and per channel with its
energies, while ours do not move when ours are refined. At the public code's second check point our number fails the
shipped check by many times its tolerance; the public code's own number there, at finer settings, lies 6.28 above its
shipped one and 0.24 from ours. Each of these is written below exactly as strongly as the measurements allow.

### 5.1 The quick example's likelihood (section 3.1's table)

What this shows: our score passes the public check by its stated tolerance, at every setting from 8 rings per panel
upward.

What this does not show: that every one of the 9,312 bin predictions equals the public code's (a single score cannot show
that; sections 5.2 and 5.3 look at the bins). Nor, from the score alone, which of the two scores, ours or the shipped one,
is nearer the exact value of the model: the shipped one is 0.093 from both of our programs. The public code run at finer
settings on 23 November 2026 answers it bin by bin, not by its score: its counts move toward ours as its mesh and its
energies are refined, and ours do not move (`CLAIMS.md`, H8 and H11); its score at this star does not approach ours
steadily (+0.093 at its shipped settings, then −0.064, −0.080, +0.029, +0.013 and +0.002 at finer ones; H10), because its
errors in energy and in its mesh partly cancel.

### 5.2 Counts, bin by bin (`quick_example_counts_check.py`)

The example ships a file of expected counts (`xpsi_good_expected_hreadable.dat`, 291 × 32 numbers). The example does not
say which run made it; X-PSI's documentation shows a synthesis at the same star with finer settings.

| | (ours / file) − 1 |
|---|---|
| total counts (577,971.4 against 577,979.8) | −1.45e-5 |
| each channel summed over the turn, channels 2 to 289 | −1.1e-5 to −1.8e-5 |
| first channel; last two channels | +1.8e-4; +4.4e-4, +1.0e-3 |
| each slice of the turn summed over channels, bright half | +1.1e-5 to +9.9e-5 |
| each slice, faint half (the spot mostly behind the star) | −6e-6 to −8.3e-4 |
| largest single difference | 1.16e-3 |

(The committed command is `quick_example_counts_check.py`, run on its own or in `check.py --slow`; the project's record
rounded these to −1.5e-5, "−1.5e-5 to −1.6e-5", +1e-5 to +1e-4, −1e-4 to −8e-4 and 1.2e-3.)

What this shows: the total and the interior channel sums agree to 2 parts in 100,000. The larger differences sit in the
first and last channels and in the faint slices.

What this does not show: why. Saying that those differences are the public code's interpolation and mesh needs the run that
made the file, its background convention and a controlled comparison. The run that made the file is still not known. The
controlled comparison of the public code's own counts with ours, at its shipped settings and finer ones, was made on
23 November 2026 (section 5.3 and `CLAIMS.md`, H8); it compares the public code as run here, not the file.

### 5.3 The public code's second shipped check (`second_check_point.py`)

`xpsi/tests/test_likelihood.py`: the same data, a different star (10 km, 10^6.8 K), no upper limit on the background,
128 energies; expected −47881.278, tolerance ±0.479.

| how | ln L | minus X-PSI's shipped value |
|---|---|---|
| ours, 48 rings per panel, 32 points, 48 energies | −47874.7600 | +6.518 |
| ours, 10 rings per panel, 16 points, 24 energies | −47874.7612 | +6.517 |
| our second program, the plain way, 48 × 192 patches | −47875.0144 | +6.264 |
| our second program, the plain way, 96 × 384 patches | −47874.8493 | +6.429 |
| X-PSI's two interpolation steps (with its 128 energies) applied to our curve (`second_check_point_emulated.py`; the record gave −47870.96 ±0.2) | −47870.95 | +10.32 |
| X-PSI itself, run here on 23 November 2026, at its shipped settings (32 cells, 64 rotation steps, 512 rays, 128 energies) | −47881.2782 | 0.000 |
| X-PSI itself, run here, at 128 cells, 256 rotation steps, 512 rays, 384 energies | −47874.9971 | +6.281 |

**Ours fails this check: +6.5 against a tolerance of 0.479.** Measured facts about this point: the star does not match the
data (the model alone gives 666,370 counts against 578,528, and exceeds the data in 238 of 291 channels), and a change of
all the model's counts by 1 part in 10,000 changes ln L by 9.4 (at the quick example's star the same change gives 0.047,
two hundred times less).

What this shows: the gap of 6.5 is what a change of all counts by 7 parts in 100,000 would give, and the tolerance 0.479
what 5 parts in a million would give. Both of our programs sit on the same side, 6.4 to 6.5 above the shipped number, and
agree with each other to 2 parts in a million in every bin; the plain one moves toward the fast one as its patches are made
smaller (−0.25, then −0.09).

One arithmetic observation (an observation, not a finding): if the public code's counts at its quick settings were higher
than ours by 6.9 parts in 100,000 in every bin, that alone would give +6.51 at this star and −0.033 at the quick example's
star (worked out with the likelihood itself; the two follow by proportion from the sensitivities of −9.44 and +0.047 per
part in 10,000 stated above). The measured gaps are +6.52 and −0.093. So one small excess of one sign matches the whole gap
here and a third of the gap there. The shipped expected-counts file, made at finer settings, is higher than
ours by 1.5 parts in 100,000 in total. This is what one would see if the public code's counts come down toward ours as its
settings are refined, and it is also what one would see by coincidence. The run at finer settings decides.

**The cause, measured on 23 November 2026: the shipped number is the public code at coarse settings for this star.** The
public code was installed and run on this laptop (`PUBLIC_CODES.md`, section 1.2; `CLAIMS.md`, H2, H3 and H9):

- at its shipped settings it gives its shipped number, −47881.278177; its count rates put through our likelihood give the
  same number to the last printed digit, so the whole gap of 6.52 lies in the counts;
- with its settings refined to 128 cells, 256 rotation steps, 512 rays and 384 energies it gives −47874.997, which is 6.28
  above its shipped number and 0.24 below ours (−47874.760); that check's tolerance is 0.479, so the public code at those
  settings would fail its own shipped check too;
- one setting at a time from its shipped ones, its number moves by +7.1 and +9.3 with its mesh (64 and 128 cells), by −2.6
  to −4.4 with its energies (192 to 512), by +1.1 with its rotation steps (128 and 256), and by +0.004 with its rays (1024
  and 2048);
- its bins move toward ours as its mesh and its energies are refined: per phase bin from 6.3e-4 at its shipped settings to
  1.2e-5 at 192 cells (with 384 rotation steps and 192 energies), per channel from 2.5e-3 to 4.8e-4 at 512 energies; ours
  do not move (three times finer settings of
  ours change ln L here by 0.000001 and no bin by more than 2.1e-7: H11).

So the arithmetic observation above is answered in part: the public code's counts do come down toward ours as its settings
are refined. But not as one excess of one sign in every bin: at its shipped settings its counts total only 3.2 parts in a
million more than ours, and the gap comes from the differences per phase bin and per channel (up to 6.3e-4 and 2.5e-3), at
a star where the score is two hundred times more sensitive to the counts than at the quick example's.

What this does not show: a value of the public code's at this star that has stopped moving. Its score is not steady across
the finer settings run (from ours: +1.85 at 64 cells, 128 rotation steps, 512 rays, 128 energies; +3.97 at 128, 256, 512,
128; +1.77 at 192, 384, 512, 192; −2.36 at 64, 128, 512, 384; −0.24 at 128, 256, 512, 384), because its errors in energy
and in its mesh partly cancel; its bins, not its score, are what move steadily toward ours. Nor does it give a time of the
public code at this star (section 6, item 13).

### 5.4 The physics

Agreement with the public code is not a result about the physics. `DERIVATION.md`, section 1, gives the four starting
points from which the programs compute light's path and says, once, what the agreement with the public code does and does
not mean for them. The result of this project is a method of calculation.

## 6. The errors found so far

The list below is the project's own list of its errors, kept whole and in the words of the one who made them: 'I' is Claude, the AI system that wrote the programs under the author's direction. Codex is a second AI system, which reviewed them.

1. Plan v0.1 said the public code makes "13 million look-ups in the hot-gas table" and would gain from the ring symmetry.
   Wrong: it already uses the ring symmetry, and the 13 million are shifted interpolations. Corrected in plan v0.2.
2. I assumed the example's channels run from 0.1 to 3.0 keV. The source shows 0.15 to 1.605 keV. Caught by the helper's
   reading before any number was made with it.
3. An observer on the spin axis got zero light (Codex found it). Corrected.
4. My fast likelihood's cut-off rule was wrong for channels with few counts (Codex's counterexample). Corrected.
5. Results v0.1 said "computes the same model", assigned the count differences to the public code, and said "no formula
   of general relativity anywhere". Each said more than was shown. Reworded.
6. My rewrite of the end-direction step first made every working-out do a slow bracketed search (the time rose from 3.2
   to 5.3 ms); found by timing, corrected the same hour.
7. My "correctness version" (`fast_star.py`) produced enormous numbers for rings seen only at their very edge. Found by a
   test against the plain way; such rings are now left out there (their true share is below 1e-7).
8. For an observer near the spin axis of a flattened star, my search for the direction in which a ring is seen edge-on
   looked outside the ring's own turn and failed (the program stopped with a message; no wrong number came out). Found by
   the geometry test against the plain way; the search now runs only between the ring's nearest and farthest directions.
9. My refusal of stars with a second path of light was incomplete: for a compact star seen from its axis with a spot near
   the far pole, the program returned one of the two paths without a word. Codex found it in round 4 with an example.
   Corrected; the example is now refused. Results v0.2 and v0.3 claimed the refusal as complete, which was wrong.
10. The compiled copy's fallback search could walk away from a root lying exactly at an end of its bracket (Codex,
    round 4). Corrected.
11. The times written at the top of the archived messages of that night ("about HH:MM") are my estimates and run up to
    1.5 hours ahead of the laptop's clock; the ledger's own times are right.
12. Before the public code was measured I expected our programs to be about 30 to 60 times faster than it in Python and
    250 to 500 times compiled, from my estimate of its time for the quick example, 0.1 to 0.2 s, scaled from its
    tutorial's time. Measured on 23 November 2026 on mains power: 49 to 55 ms at its shipped settings, so 17 times in
    Python and 116 to 131 times compiled. My estimate of its time was about 2 to 4 times too long, and the ratios I expected
    about 2 to 4 times too high.
13. My stopwatch for the public code (`time_public_code.py`) times it at stars near the check's star, each number moved by
    up to 1 part in 1000. At the second check that star sits on the public code's lower limit of the radius, so about half
    of the nearby stars were refused at once (3 of 4, and 1 of 2, at the two settings that timed more than one star), and
    their times are not times of an evaluation. Found when it was first run. Times and ratios are now given for the quick
    example only (`compare_with_public_code.py` leaves refused stars out); the values and the bins at the check's star
    itself are not affected.
