# Derivation: what the programs compute, and from what

This page is for a reader who doubts and wants to check: a physicist or a programmer, not an insider. It derives, step by step and from named starting points, what the two Python programs and the compiled copy compute. Each statement is given first in plain words and then exactly. Each formula names the function that computes it (file and function) and the check that tests it numerically. The quick checks run with `python -I -B -X utf8 2-python-solution/check.py`, the slow ones with the same command and `--slow`; the check names quoted on this page, in double quotes, are the names that command prints. The claims table `CLAIMS.md` in this folder lists the numbers stated on this page with where each was measured.

Files named on this page, all under `2-python-solution/code/` unless another folder is given: `bend.py` (the path integrals), `fast_kernel.py` (the fast program), `plain_star.py` (the slow program, independent in method), `public_likelihood.py` (the public code's likelihood, written anew), `quick_example.py` (the public code's shipped example worked with our programs), `consts.py`; `2-python-solution/tools/` (the slow checks, each its own script); `3-compiled-solution/nicer_fast.c` (the compiled copy of the fast program, with its own check command `3-compiled-solution/check.py`).

Notation used throughout: `r` is the radius (the circumference of the circle at that radius, over 2 pi); `L` is the length in the clock rate; `f(r) = sqrt(1 - L/r)`; `alpha` is the angle between a ray and the outward radius as the frame at rest there reads it; `x = cos(alpha)` at the surface; `psi` is an angle swept about the centre; `t` is far time (the time of clocks far from the star; the observer's time); `c` is the speed of light. The star: `Req` its equatorial radius, `R` (or `R(theta)`) its radius at colatitude `theta`, `Omega = 2 pi nu` its spin rate in far time, `D` the distance to the observer, `i` the observer's colatitude (the public example gives `cos_inclination`). The public example's parameter names are kept where they are meant: mass, radius, distance, cos_inclination, phase_shift, super_colatitude, super_radius, super_temperature.

## 1. Starting points

Labels (the project's rule): DERIVED means derived in the author's second paper ("Reconciling the Schwarzschild Solution with Conservation of Mass-Energy", draft v0.46); this page cites that paper and does not re-derive it. IMPORT means a law not in the author's papers, named here as a premise. INPUT OF THE PUBLIC EXAMPLE means a definition the example makes of its star, used because the example defines its star so.

**P1, the clock rate.** Plain: a clock held at rest at radius r outside a round body that does not spin runs slow against clocks far away, by the factor f(r). Exact: `f(r) = sqrt(1 - L/r)`; while far clocks advance `dt`, a clock at rest at r advances `dt_loc = f(r) dt`. Label: DERIVED in the author's second paper, section 3.

**P2, the rulers.** Plain: a ruler laid along the radius measures more than the change of r; a ruler laid sideways measures r times the angle. Exact: a radial step `dr` has ruler length `dr / f(r)`; a sideways step `dpsi` has ruler length `r dpsi`, r being the circumference over 2 pi. Label: DERIVED there (the ruler rule).

**P3, light's speed.** Plain: in every frame at rest, light moves at c. Exact: the paper's H2, special relativity holds in each frame at rest; light moves at c in each such frame. Label: DERIVED there. Where this page uses more of special relativity than the speed of light (the Doppler factor, aberration, time dilation, length contraction), it uses the Lorentz transformation in the frame at rest at the emission point, as a consequence of H2, with its steps shown (section 3.2), and nothing else.

**P4, the constant of the path.** Plain: along a ray, the radius times the sine of the angle to the radius, divided by the clock rate, does not change. Exact: `b = r sin(alpha) / f(r)` has the same value at every point of a ray. Label: IMPORT, accepted by the author as a premise on 23 November 2026. This law (the angular momentum law for light) is not derived in the author's two papers. Every step on this page that uses it says so.

**Arrival.** Plain: a photon that leaves the surface in an open direction reaches the far observer; nothing in between absorbs it. Label: IMPORT. "Open" means above the surface's tangent plane (`mu > 0` in section 3.2), along the one path of light the program finds (section 3.4 says what the program refuses).

**Inputs of the public example, flagged.**

- The star's flattened shape: `R(theta) = Req (1 + o2 cos^2 theta)`, `o2 = (Omega^2 Req^3 / (G M)) (-0.788 + 1.030 G M / (Req c^2))`. The public example defines its star this way, and that formula was made to match star models of general relativity. It is used only because the example defines its star so. It is the line computing `o2` in `fast_kernel.spot_counts` and in `plain_star.light_curve` (and in `nicer_fast.c` `spot_counts`); with `oblate=False` the star is round, as the test article's round-star curves define it.
- The length L from the star's mass: the programs set `L = 2 G M / c^2` (the line computing `L` in `fast_kernel.spot_counts` and `plain_star.light_curve`), with the public code's own value of G M_sun / c^2 for the example (`quick_example.star_args`) and the observed G M_sun of `consts.py` for the test article's stars. The relation itself is DERIVED in the author's second paper (section 3, step 4: the count inside a closed surface fixes the clock rate, `f^2 = 1 - 2 G M / (c^2 r)`, so `L = 2 G M / c^2`). What the public example and the test article supply are the values put into it: the star's mass and the constant.
- The spot: a cap of angular radius super_radius about the colatitude super_colatitude, radiating as a black body at the temperature super_temperature, the same in every direction in its own rest frame: `N'(E') = BB_CONST E'^2 / (exp(E'/kT) - 1)` photons per unit area, time, solid angle and energy, times the cosine of the angle to the surface's normal (`fast_kernel.bb_photons`; for the example the public code's own black-body constant, `X_BB` in `quick_example.py`). Four of the test article's curves add a beaming function `g` of that cosine (`beaming` in `fast_kernel.spot_counts`). These are definitions of the example's and the test problems' star, not physics of ours.
- The paths outside the flattened, spinning star. P1 to P4 are stated for a round body that does not spin. The programs use them outside the flattened, spinning star exactly as for a round body with the same L, from each emission point's own radius `R(theta)`: the flattening and the spin enter only at the surface (where the light leaves, how the surface is tilted, how fast it moves). That this is right for a spinning, flattened body is not derived on this page: step not yet written out.

For a round body that does not spin, P1 to P4 give the same light paths the public code uses, so this project cannot tell the two apart. Agreement with the public code is never a proof of the physics and never a result about the physics. The result of this project is a method of calculation.

## 2. From P1 to P4 to the two path integrals

### 2.1 In the frame at rest

Plain: follow a ray for a moment in the frame at rest where it is. There it moves at c, so its radial and sideways advances are the two legs of a right triangle with hypotenuse c times the frame's time. The rulers turn those legs into changes of r and of the angle; the clock turns the frame's time into far time.

Exact. At radius r the frame at rest reads the ray's direction as the angle `alpha` to the outward radius. In that frame light moves at c (P3), so in the frame's own time `dt_loc` the ray advances `c dt_loc`, of which `c cos(alpha) dt_loc` is radial and `c sin(alpha) dt_loc` sideways, both ruler lengths. By P2 the ruler length of a step `dr` is `dr / f` and that of a turn `dpsi` is `r dpsi`:

    dr / f = c cos(alpha) dt_loc,        r dpsi = c sin(alpha) dt_loc.

By P1, `dt_loc = f dt`. Dividing the second equation by the first:

    dpsi / dr = sin(alpha) / (r f cos(alpha)).

By P4 (IMPORT), `sin(alpha) = b f / r` all along the ray, so `cos(alpha) = sqrt(1 - b^2 f^2 / r^2)`. Write

    Q = 1 - b^2 f^2 / r^2 = cos^2(alpha).

Then `dpsi / dr = (b f / r) / (r f sqrt(Q))`, that is

    dpsi / dr = b / (r^2 sqrt(Q)).

From the first equation with `dt_loc = f dt`: `dr = f c cos(alpha) dt_loc = f^2 c cos(alpha) dt`, so

    c dt / dr = 1 / (f^2 sqrt(Q)).

These two lines are the "CONSEQUENCES" in the docstring of `bend.py`. One geometric step is used alongside them: the ray stays in the plane that holds the centre, the emission point and its starting direction. Outside a round body that does not spin, that plane is a plane of mirror symmetry of the body and of the ray's start; a ray that left the plane would have a mirror image that is a second ray with the same start, and the equations above give one ray per start. So `psi` is the angle swept in that plane.

### 2.2 Units, and the two integrals

Plain: measure lengths in the star's radius at the emission point, and the delay in the time light needs for that radius. Then the swept angle and the delay are two one-dimensional integrals over the same square root, with two numbers as inputs: how compact the star is, and the direction the light leaves in.

Exact. Let `y = R / r` (1 at the surface, 0 far away), `u = L / R`, so `f^2 = 1 - L / r = 1 - u y`. Let `x = cos(alpha)` at the surface. P4 (IMPORT) at the surface gives `b = R sin(alpha) / f(R)`, so

    kb^2 = (b / R)^2 = (1 - x^2) / (1 - u).

Then `b^2 f^2 / r^2 = kb^2 R^2 (1 - u y) y^2 / R^2`, and

    Q = P(y) = 1 - kb^2 y^2 (1 - u y).

With `r = R / y`: `dr = -R dy / y^2` and `1 / r^2 = y^2 / R^2`, so `b dr / r^2 = -(b / R) dy`; r from R to far away is y from 1 to 0, and the minus sign turns the limits round. The swept angle from the surface to far away:

    Theta(u, x) = INT_R^inf (dpsi/dr) dr = INT_R^inf b dr / (r^2 sqrt(Q)) = kb INT_0^1 dy / sqrt(P).

The delay. The far time a ray needs from the surface to a far radius `r_far` is `c t_ray = INT_R^r_far dr / (f^2 sqrt(Q))`. A ray that leaves the same point straight outward has `alpha = 0`, `b = 0`, `Q = 1`, and needs `c t_rad = INT_R^r_far dr / f^2`. Each grows without bound as `r_far` grows (both integrands tend to 1), but their difference settles, because `1 / sqrt(Q) - 1` falls off as `1 / r^2`:

    c (t_ray - t_rad) = INT_R^inf (1 / sqrt(Q) - 1) dr / f^2.

This is why the delay's integrand is `1 / (f^2 sqrt(Q))`, the ray's own `c dt/dr`, minus the radial ray's `1 / f^2`. In y, `dr / f^2 = -R dy / (y^2 (1 - u y))`, so in units of `R / c`

    tau = INT_0^1 (1 / sqrt(P) - 1) dy / (y^2 (1 - u y)).

Far away (`y` near 0) `P` is near 1 and `1 / sqrt(P) - 1` is the difference of two numbers near 1, which loses digits. Multiply above and below by `1 + sqrt(P)`: `1 / sqrt(P) - 1 = (1 - sqrt(P)) / sqrt(P) = (1 - P) / (sqrt(P) (1 + sqrt(P)))`, and `1 - P = kb^2 y^2 (1 - u y)` cancels the denominator `y^2 (1 - u y)` exactly:

    tau(u, x) = kb^2 INT_0^1 dy / (sqrt(P) (1 + sqrt(P))).

Near the surface (`y` near 1) `P` falls to `x^2`, which for a ray leaving nearly along the surface is small and would again be a difference of numbers near 1. With `w = 1 - y`: `y^2 = 1 - 2w + w^2`, `y^3 = 1 - 3w + 3w^2 - w^3`, and `P - x^2 = kb^2 (1 - u) - kb^2 y^2 (1 - u y) = kb^2 [(1 - u) - y^2 + u y^3] = kb^2 [2w - w^2 - 3uw + 3uw^2 - uw^3]`, so exactly

    P = x^2 + kb^2 w ((2 - 3u) - (1 - 3u) w - u w^2),

a sum of `x^2` and a term that is zero at the surface and positive just inside the range of integration (for `u < 2/3`). The programs integrate in `v = sqrt(w) = sqrt(1 - R/r)`, `dy = -2 v dv`, and where x is small they use `v = (x / sqrt(a)) sinh(s)`, `a = kb^2 (2 - 3u)`, which keeps the integrand smooth (`bend.theta_tau` for `x < 0.25`; `fast_kernel.bend` for every x). These are substitutions inside the integral, not steps of physics.

Functions: `bend.theta_tau` (outward rays, Gauss-Legendre with 160 points), `fast_kernel.bend` (any direction, see 2.4). In the compiled copy the integrals are not evaluated at run time; see 2.5.

### 2.3 The common reference: a radial ray from the equator

Plain: the delays of all rings are referred to one ray, the one that leaves the equator straight outward; a ring at another latitude sits at another radius, and the radial ray from there needs a different time to reach the equatorial radius.

Exact. The programs use `dly = Omega (tau R / c + trad)` with

    trad = (Req - R + L ln((Req - L) / (R - L))) / c

(the lines computing `trad` in `fast_kernel.spot_counts` and `plain_star.light_curve`). For the radial ray `c dt/dr = 1 / f^2 = r / (r - L) = 1 + L / (r - L)`, and `INT_R^Req (1 + L / (r - L)) dr = (Req - R) + L ln((Req - L) / (R - L))`. The far time of the radial ray from Req to the observer is the same for every photon of the star and is dropped: it shifts all phases alike, and the public example fixes its phase convention with phase_shift (section 4).

### 2.4 Rays that start inward (x < 0)

Plain: on a flattened star the surface is tilted against the sphere, so a direction that points a little toward the centre can still lie above the surface. Such a ray first goes in, turns where it runs sideways, and comes back out through the surface's radius; from there it is an ordinary outward ray.

Exact. The equations of 2.1 depend on r only, so between the turning radius `r_t` (where `Q = 0`, that is, the frame there reads `alpha = 90 degrees`) and R the inward and the outward halves sweep the same angle and take the same far time. Hence, for `x < 0`,

    Theta(u, x) = Theta(u, |x|) + 2 (angle swept between r_t and R),
    tau(u, x)   = tau(u, |x|)   + 2 (far time between r_t and R, in units of R / c);

on the extra stretch no radial ray is subtracted, since it is extra path. With `s = R / r - 1 >= 0` on the stretch, `y = 1 + s`, and `(1 + s)^2 (1 - u - u s) = (1 - u) + (2 - 3u) s + (1 - 3u) s^2 - u s^3`, so with `c1 = 2 - 3u`, `c2 = 1 - 3u`, `c3 = -u`:

    P = 1 - kb^2 (1 + s)^2 (1 - u (1 + s)) = x^2 - kb^2 s (c1 + c2 s + c3 s^2)

(using `kb^2 (1 - u) = 1 - x^2`). The turning point `s_t` is the root of `kb^2 s_t (c1 + c2 s_t + c3 s_t^2) = x^2`, found by Newton's steps. Subtracting the two cubics,

    P = kb^2 (s_t - s) g(s),    g(s) = c1 + c2 (s_t + s) + c3 (s_t^2 + s_t s + s^2),

which has no cancellation. With `dr = -R ds / (1 + s)^2` and `f^2 = 1 - u (1 + s)`, the extra angle is `2 kb INT_0^s_t ds / sqrt(P)` and the extra far time `2 INT_0^s_t ds / ((1 + s)^2 (1 - u (1 + s)) sqrt(P))`; the substitution `s = s_t (1 - v^2)` removes the square-root end. Functions: `bend.theta_tau_any`, `fast_kernel.bend`.

### 2.5 The approximation of the path integrals

Plain: the swept angle and the delay depend on two numbers only, u and x. Both are smooth after a square root is taken out, so one table of polynomial coefficients, built once from the integrals themselves, serves every star in the range the programs cover.

Exact. Near `x = 1` (the radial ray) `kb^2 = (1 - x)(1 + x) / (1 - u)`, so `Theta` is proportional to `z = sqrt(1 - x)` and `tau` to `z^2`; `Theta / z` and `tau / z^2` are smooth in `z` and in `u`. The approximation of the path integrals is a two-variable Chebyshev expansion of these two functions, with `FIT_NU = 32` points in `u` on `FIT_U = (0.15, 0.55)` and `FIT_NZ = 72` points in `z` from 0 to `sqrt(1 - FIT_XMIN)`, `FIT_XMIN = -0.12`, built once from `fast_kernel.bend` with 96 Gauss points (`fast_kernel._path_fit`). For each ring the two-variable coefficients are contracted in `u` to a one-variable series in `z` (`_ring_series`), evaluated at the ring's directions (`_series_eval`); in the compiled copy `swept` and `swept_delay` in `nicer_fast.c` evaluate those series. The program uses the approximation when the star's `L / R` lies inside `FIT_U` at every latitude and every ring's tilt has `|sin|` below `FIT_TILT = 0.08`; otherwise it uses the integrals (`fit=True` in `fast_kernel.spot_counts` chooses between them). Measured against `fast_kernel.bend` with 96 points on 6,000 random directions (`2-python-solution/tools/path_approximation_check.py`, in `check.py --slow`: "the approximation of the path integrals against the integrals themselves on 6000 random directions"): largest difference 8.5e-13 rad in the swept angle and 1.6e-12 R/c in the delay (the comment in `fast_kernel.py` gives 9e-13 rad and 1.3e-12 R/c from another sample of directions).

### 2.6 An identity of the two integrals, and the three checks of light's path

Plain: moving the light's source a little sideways changes the arrival time by the part of that move that lies along the ray, counted in far time. The two integrals must obey this; it is a check that depends on neither program.

Exact. At fixed u, both integrals depend on x through `kb` alone. Differentiate under the integral sign (for `x > 0`, `P >= x^2 > 0` on the whole range, so this is allowed). Since `dP / dkb = -2 kb y^2 (1 - u y)`:

    d(P^(-1/2)) / dkb = kb y^2 (1 - u y) / P^(3/2).

For the swept angle, `Theta = kb INT dy / sqrt(P)`:

    dTheta / dkb = INT dy / sqrt(P) + kb^2 INT y^2 (1 - u y) dy / P^(3/2) = INT [P + kb^2 y^2 (1 - u y)] dy / P^(3/2) = INT_0^1 dy / P^(3/2),

because `P + kb^2 y^2 (1 - u y) = 1`. For the delay, in the form `tau = INT (1 / sqrt(P) - 1) dy / (y^2 (1 - u y))` of 2.2, only `P` depends on `kb`:

    dtau / dkb = INT kb y^2 (1 - u y) dy / (P^(3/2) y^2 (1 - u y)) = kb INT_0^1 dy / P^(3/2) = kb dTheta / dkb.

Hence `dtau / dTheta = kb = b / R`; in units of R for b and of R/c for the delay T, `c dT / dTheta = b`. In words: a sideways move of the emission point by ruler length `R dTheta` changes the far arrival time by `sin(alpha) R dTheta / (c f)`, the component of the move along the ray, divided by c and by the clock rate (P1, P2, P4 IMPORT). The same step shows `dTheta / dkb > 0`: the swept angle grows as x falls, which section 3 uses.

The three checks in `2-python-solution/check.py`, on the test article's star of 1.4 Suns and 12 km:

- "swept angle against the article's table (its printed precision)": `fast_kernel.bend` with 96 points compared with the article's table of swept angle against leaving angle, for leaving angles below 1.5607 rad; limit 2e-5, relative.
- "delay behind a radial ray against the article's table (its printed precision)": `bend.theta_tau` times `R / c` compared with the article's table of delays; limit 2e-5, relative.
- "the delay's rate of change with the swept angle against b (an identity of the two integrals)": `(tau(x + h) - tau(x - h)) / (Theta(x + h) - Theta(x - h))` compared with `sqrt((1 - x^2) / (1 - u))` at 49 points x from 0.02 to 0.98, `h = 1e-5`; limit 1e-7. This check needs no table and no second program.

## 3. The two exact results the fast program rests on

These two are identities, not approximations; both were derived independently and checked line by line in the programs. What is approximate is the arithmetic around them: a finite number of rings, of points and of energies (3.4).

### 3.0 The ring and its turn

Plain: cut the spot into rings of latitude. A point of a ring goes round once per turn; its angle from the observer's direction swings between a nearest and a farthest value. Light that leaves it in the right direction reaches the observer.

Exact. A ring at colatitude `theta`; a point of it at turning angle `phi`, measured from the meridian that holds the observer's direction (`phi = 0` when the point is nearest to it). The angle `psi` between the point's radius and the observer's direction follows the spherical cosine rule:

    cos(psi) = A + B cos(phi),    A = cos(i) cos(theta),    B = sin(i) sin(theta).

The light that reaches the observer leaves in the direction x with `Theta(u, x) = psi`. On the half of the turn before closest approach (`s = -1`) and on the half after it (`s = +1`), `sin(phi)` has one sign and x is a monotone function of `phi`, because `psi` is monotone in `phi` on each half and `Theta` is monotone in x (2.6). The ends of a ring's range of x: `x_hi`, the direction for the nearest point (`psi = |theta - i|`; `x_hi = 1` when that angle is 0), and `x_lo`, either the direction for the farthest point (`psi = theta + i`, or `2 pi - theta - i`) when the whole turn is in view, or the edge-on direction where `mu = 0` (3.2) when part of the turn is hidden; for a round star the edge-on direction is `x = 0` exactly. They are found by secant steps on `Theta(u, x) = psi` with the residual checked and bisection where a step fails (`solve_rows` inside `fast_kernel.spot_counts`). The spot's stretch of the ring is the longitude range `[-dl, dl]` about the spot's centre, `cos(dl) = (cos(super_radius) - cos(theta) cos(super_colatitude)) / (sin(theta) sin(super_colatitude))` (`fast_kernel.spot_rings`). The rings are Gauss points in latitude on panels; the panels are cut at the spot's edges, at the latitude where the spot starts to cover whole rings (a pole inside it), and at the latitudes where a ring's visibility changes kind (part of the turn becomes hidden; the whole turn comes back into view; the ring is never seen beyond), because a ring's contribution has a corner there with square-root behaviour on one side, and a panel end with that behaviour gets a squared variable (`fast_kernel.visibility_breaks`, `spot_rings`, `_panel`). Rings whose turn does not change `psi` (`B` near 0: the observer on the spin axis, or the ring at a pole) give a steady flux and are added separately, with a lensing factor by a centred difference (the `aligned` branch of `spot_counts`); everything below is for `B > 0`.

### 3.1 Result a: counting the photons at emission removes the lensing factor

Plain: the usual way takes each moment of the turn, finds the direction that reaches the observer (an inversion), and multiplies by how much a bundle of directions is spread by the bending (the lensing factor, a derivative). Counting instead by direction, the spreading factor cancels against the change of variable, exactly: for each direction, the moment of the turn follows by the cosine rule, forward.

Exact. Along the turn, `cos(psi) = A + B cos(phi)` gives `d(cos psi) = -B sin(phi) dphi`. Also `cos(psi) = cos(Theta(u, x))`, so `d(cos psi) = -sin(psi) Theta'(x) dx`, with `Theta' = dTheta/dx`. The lensing factor is

    J = |d cos(alpha) / d cos(psi)| = |dx / d(cos psi)| = 1 / (sin(psi) |Theta'(x)|).

Equating the two expressions for `d(cos psi)`: `|dx| = B |sin(phi)| |dphi| / (sin(psi) |Theta'(x)|)`, so `|dphi| / (sin(psi) |Theta'(x)|) = |dx| / (B |sin(phi)|)`, and on each half of the turn (where `sin(phi)` keeps one sign)

    J |dphi| = |dx| / (B sin(phi)).

So an integral over the turn of (anything) times `J dphi` is an integral over x of (anything) `/ (B sin(phi))` dx, where for each x the program computes `psi = Theta(u, x)` forward, then `cos(phi) = (cos(psi) - A) / B` and `sin(phi)`. No inversion of `Theta`, no derivative of it. P4 (IMPORT) enters through `Theta`.

Function: `fast_kernel.spot_counts`. The directions are Chebyshev points `x = xm + xh cos(t)`, `xm = (x_hi + x_lo) / 2`, `xh = (x_hi - x_lo) / 2`, `t` in `[0, pi]` (`_cheb`), so `|dx| = xh sin(t) dt`, and the weight `wg` carries `xh sin(t) / (B sin(phi))`. Checks:

- `2-python-solution/tools/count_at_emission_check.py`: for one small surface element on the test article's star, the usual way (a table of directions, its inversion, the lensing factor by differences of the table, 400,000 emission phases) compared with counting at emission, through the Fourier coefficients of the element's light over one turn, at 8, 16, 32 and 64 nodes.
- "fast kernel against the plain way, bin by bin, of the mean bin (plain way at coarse patches)" in `check.py`: `spot_counts` compared with `plain_star.light_curve` on the example's star at four energies; limit 2e-4 of the mean bin. The plain program (`plain_star.py`) keeps the usual way: a mesh of patches over the spot, a fine grid of emission phases, for each patch its own table of direction against angle, the lensing factor by differences of that table, the flux resampled at arrival. It shares with the fast program only the path integrals and the inputs.
- The slow checks: "round stars, 36 geometries against the plain way: largest difference of the largest bin" (limit 4e-4; `kernel_against_plain_way.py`; measured: 34 geometries within 2e-4 of the largest bin, most near 1e-5, two at 2.4e-4 and 3.2e-4 where only a sliver of the spot is in view; that those two fall to 1.1e-4 and 2.1e-4 at four times as many patches of the plain way is the project's record, with no committed command), and "quick example's star: fast kernel (48 rings per panel, 32 points) against the plain way at 96 x 384 patches, bin by bin, of the mean bin" (limit 1e-5; `likelihood_from_plain_star.py`; measured: 3.7e-6, 4 parts in a million, in every bin; on the second check star 2.1e-6). On its own command, `kernel_against_plain_star_many.py` (about fifteen minutes): 24 flattened-star geometries at 200 to 716 Hz, 22 worked out, all within 1.9e-5 of the largest bin, 1 refused by the fast program (a second path of light is possible), 1 by the plain way.
- Against the agreed reference curves of the test article: "fast kernel against the agreed reference curve OS1e, of its median (the article's target 1e-3)" (quick), and the slow "the twelve agreed reference curves with extended spots, 24 rings per panel and 24 points, of the median (the article's target)" (limit 1e-3; measured 1.2e-5 to 2.8e-4), "the same twelve curves at the working setting (10 rings per panel, 16 points), of the median" (limit 1e-3), "the four curves with point-like spots (radius 0.01), of the median: above the article's target at eclipse edges" (limit 1e-2; measured 2.4e-3 to 7.7e-3, at the edges of eclipses) (`kernel_against_reference_curves.py`).

### 3.2 The count per emitted photon, and why the two programs agree

Plain: a patch of the hot surface emits, in its own frame, a black body's photons. The patch moves, so the frame at rest at that point sees its photons shifted in energy and direction (Doppler); the clock rate shifts the energy once more on the way out; the bending spreads the directions (the lensing factor). The fast program counts the photons as they leave, per moment of the turn; the slow program counts them as they arrive, per second of arrival. The two counts differ by the ratio of arrival time to emission time, and that ratio is exactly the Doppler factor's companion, through the identity of 2.6.

Exact. Let S be the frame at rest at the emission point, S' the frame of the surface element, which moves relative to S at speed `v = beta c` along the ring (the azimuthal direction, tangent to the surface and perpendicular to its normal).

The speed. In far time `dt` the element turns by `dphi = Omega dt`; the ruler length of that arc is `R sin(theta) dphi` (P2, sideways rulers on the circle of latitude of circumference `2 pi R sin(theta)`); the frame's own time is `dt_loc = f dt` (P1). So

    beta = Omega R sin(theta) / (c f),    gamma = 1 / sqrt(1 - beta^2)

(the lines computing `beta` and `gam` in `fast_kernel.spot_counts` and `plain_star.light_curve`).

The Lorentz transformation (H2). For a photon of energy `E` and momentum `p = E / c` (light, P3) at angle `xi` to the velocity in S, special relativity in S gives in S': `E' = gamma (E - v p_par)`, `p'_par = gamma (p_par - v E / c^2)`, `p'_perp = p_perp`, with `p_par = (E / c) cos(xi)`; a clock moving with S' advances `dt'` while S's clocks advance `gamma dt'`; a length along v that measures `l'` in S' measures `l' / gamma` in S. Five consequences:

1. Energy: `E' = gamma E (1 - beta cos(xi))`, so `E = delta E'` with
   `delta = 1 / (gamma (1 - beta cos(xi)))`.
2. Direction: `cos(xi') = p'_par c / E' = (cos(xi) - beta) / (1 - beta cos(xi))`. Differentiating, `d cos(xi') = [(1 - beta cos xi) + beta (cos xi - beta)] d cos(xi) / (1 - beta cos xi)^2 = (1 - beta^2) d cos(xi) / (1 - beta cos xi)^2 = delta^2 d cos(xi)`; the azimuth about v is unchanged, so a solid angle `dOmega_S` in S is `dOmega' = delta^2 dOmega_S` in S'.
3. The normal `n` is perpendicular to v, so `p'_n = p_n`: `E' mu' = E mu`, where `mu`, `mu'` are the cosines of the angle to the normal in S and S'. Hence `mu' = delta mu` (the argument of the beaming function).
4. Area: the element's area `dS'` in its own frame is `gamma` times its area `dS_S` in S (contraction along v).
5. Time: the element's own time `dt'` is `dt_S / gamma`.

The angle `xi` in S. The ray's direction has outward component `cos(alpha)` and sideways component `sin(alpha)` along the great circle from the point toward the observer's direction. The cosine of the angle between the ring's tangent and that great-circle direction is `-sin(i) sin(phi) / sin(psi)` (from `d(cos psi)/dphi = -B sin(phi)` and the arc `R sin(theta) dphi`), so

    cos(xi) = -sin(alpha) sin(i) sin(phi) / sin(psi)

(`cxi` in `plain_star.light_curve`; `k1 = |cos xi|` in `fast_kernel.spot_counts`, with the sign set by the half-turn: approaching before closest approach).

The normal and `mu`. Along the meridian, a step `dtheta` moves `R dtheta` sideways and `dR = R'(theta) dtheta` in radius, of ruler length `dR / f` (P2), so the surface's slope against the sphere is `h = R'(theta) / (R f)` with `R' = -2 Req o2 cos(theta) sin(theta)`; the element's area is `R^2 sqrt(1 + h^2) dOmega_patch` (`dOmega_patch` the patch's solid angle from the centre), and the normal has components `cg = 1 / sqrt(1 + h^2)` along the radius and `sg = h / sqrt(1 + h^2)` toward the pole. The cosine of the angle at the point between the direction toward the pole and the direction toward the observer is `(cos(i) - cos(theta) cos(psi)) / (sin(theta) sin(psi))` (the spherical cosine rule in the triangle pole, point, observer's direction). So

    mu = x cg + sin(alpha) sg (cos(i) - cos(theta) cos(psi)) / (sin(theta) sin(psi)),

and the light leaves only where `mu > 0` (the lines computing `h`, `cg`, `sg`, `mu` in both programs). For a round star `h = 0` and `mu = x`.

The energy shift to far away. A clock at rest at R advances `f dt` while far clocks advance `dt` (P1); the crests of a wave are counted the same by both, so a frequency `nu_S` in S is `f nu_S` far away: `E_obs = f E_S = f delta E'`, and the emitted energy that arrives as `E_obs` is

    E' = E_obs / (f delta)

(the argument of `bb_photons` in both programs).

The arrival interval. Photons that leave the moving element during far time `dt_em` arrive during `dt_arr = dt_em + d(delay)`. The delay is `(R / c) tau + trad` with `trad` fixed on the ring, so `d(delay)/dphi = (R / c) (dtau / dTheta) (dpsi / dphi) = (R / c) kb (B sin(phi) / sin(psi))` by 2.6, and with `kb = sin(alpha) / f` (from `kb^2 = (1 - x^2) / (1 - u)`), `B = sin(i) sin(theta)` and `dphi = Omega dt_em`:

    d(delay) = (Omega R sin(theta) / (c f)) (sin(alpha) sin(i) sin(phi) / sin(psi)) dt_em = -beta cos(xi) dt_em,
    dt_arr = (1 - beta cos(xi)) dt_em = dt_em / (gamma delta).

The spreading of directions. Rays that leave the point at the same `alpha` in different azimuths about the radius each stay in their own plane (2.1) and reach the far sphere at the same angle `psi = Theta(u, x)` from the point's radius. So a solid angle `dOmega_S` of leaving directions in S becomes `dOmega_far = dOmega_S |d cos(psi) / d cos(alpha)| = dOmega_S / J` far away, and the observer's collecting area `dA_obs` at distance D (far: `f(D)` taken as 1) sees `dOmega_far = dA_obs / D^2`. All photons that leave arrive (IMPORT).

Assembling the count. In S', during `dt'`, the element emits into `dOmega'` about a direction at cosine `mu'` to the normal, in the energy range `dE'`:

    dN = N'(E') g(mu') mu' dE' dS' dOmega' dt'

(the black-body spot, an input of the example; `g = 1` without beaming). Replace each factor by its value in S, using the five consequences: `N'(E_S / delta) g(delta mu) (delta mu) (dE_S / delta) (gamma dS_S) (delta^2 dOmega_S) (dt_S / gamma)`, that is

    dN = delta^2 mu N'(E_S / delta) g(delta mu) dE_S dS_S dOmega_S dt_S.

Now to the far observer: `dE_S = dE_obs / f`, `dt_S = f dt_em` (P1), `dOmega_S = J dA_obs / D^2`; the two factors f cancel:

    dN = (dS_S / D^2) J delta^2 mu N'(E_obs / (f delta)) g(delta mu) dE_obs dA_obs dt_em.

Per unit far time of emission this is the count the fast program adds up: with `dS_S = R^2 sqrt(1 + h^2) sin(theta) dtheta dlambda`, `J dphi = dx / (B sin(phi))` (3.1) and `dphi = Omega dt_em`, and with the arrival bins measured in turns (one turn is `2 pi` of arrival angle; the Omega cancels), the weight per unit `t` (the Chebyshev variable), per unit latitude and longitude of the stretch, is

    w[r, s, j, E] = (xh sin(t) / (B sin(phi))) (delta^2 / gamma) mu N(E / (f delta); delta mu)

times `gamma R^2 sqrt(1 + h^2) sin(theta)` (the area with its gamma), divided by `D^2` and by `2 pi`: `wg`, `Np` and the lines `area = ...` and `Mk *= area * wth / (D_cm * D_cm * two_pi)` in `fast_kernel.spot_counts`; `spot_counts` in `nicer_fast.c`. The `gamma` of the area and the `1 / gamma` of the weight cancel; they are kept apart because the test article's round-star curves are defined without the gamma of the area (`gamma_area=False`).

Per unit arrival time, divide by `dt_arr / dt_em = 1 / (gamma delta)`:

    dF = (dS' / D^2) J delta^3 mu N'(E / (f delta)),    dS' = gamma R^2 sqrt(1 + h^2) dOmega_patch,

the plain program's flux per patch (`geom` and `F` in `plain_star.light_curve`). The two agree because `delta^3 dt_arr = delta^3 dt_em / (gamma delta) = (delta^2 / gamma) dt_em`: the identity `c dT / dTheta = b` of 2.6 is what the agreement rests on. Where P4 (IMPORT) enters: in `Theta`, in `J`, and in `dt_arr` through the identity.

Checks: the ones of 3.1 compare the whole count, fast against plain, bin by bin; the two programs share the path integrals, the shape, the tilt, the area, the energy shift, the Doppler factor and the factors above, so their agreement checks the adding-up and the identity, not those shared formulas. The agreed reference curves of the test article (3.1, last item) check the shared formulas, including the beaming argument `delta mu` (the curves OS1g to OS1j).

### 3.3 Result b: a whole stretch of the spot goes into a phase bin exactly

Plain: every point of a ring's stretch runs through the same moments of the turn, only earlier or later by its longitude. So the stretch's share of an arrival bin is the share of one point, averaged over a window of longitudes; and an average over a window of a count is a second difference of the count added up once more. No patch across the stretch, no curve to interpolate.

Exact. Let `w(phi) dphi` be the count of 3.2 for the point of the ring at the stretch's middle, per unit longitude of the stretch, and `a(phi) = s phi + Omega delay(phi)` its arrival angle (in the program `a = dly -/+ phi`, with `dly = Omega (tau R / c + trad)`; `phase0` fixes where the spot's centre is nearest). A point of the stretch at longitude `lambda` ahead reaches each turning angle earlier by `lambda / Omega`, so its photons emitted at `phi` arrive at `a(phi) - lambda`. The count of the stretch that lands in the arrival bin `[b0, b1]` is

    INT_{-dl}^{dl} dlambda INT dphi w(phi) [b0 <= a(phi) - lambda < b1].

Define `N(z) = INT dphi w(phi) [a(phi) < z]`, the count of emitted photons whose arrival angle is below z, and

    K(z) = INT_{-inf}^{z} N(z') dz' = INT dphi w(phi) max(z - a(phi), 0),

the count of arrived photons, added up once more: the sum over emitted photons with `a < z` of `(z - a)`. The inner integral over `phi` is `N(b1 + lambda) - N(b0 + lambda)`, and `INT_{-dl}^{dl} N(b + lambda) dlambda = K(b + dl) - K(b - dl)`, so the share is

    K(b1 + dl) - K(b0 + dl) - K(b1 - dl) + K(b0 - dl).

This holds for any functions `w` and `a`. In the program `a(phi)` is monotone on each half-turn (`da / dphi = 1 - beta cos(xi) > 0`, 3.2, and `phi` is monotone in x), so `{a < z}` is one interval of the Chebyshev variable t, from the half-turn's start to the point `t*` where `a = z`. The integrand `w (z - a)` is represented by its Chebyshev interpolant through the n nodes, `sum_j w_j (z - a_j) l_j(t)`, with `l_j` the cardinal function of node j, and

    K(z) = sum_j w_j (z - a_j) p0_j(t*),    p0_j(t*) = INT_0^{t*} l_j(t) dt,

where the `p0_j` are polynomials with known coefficients (the matrix `G` of Chebyshev primitives; `q0 = p0_j(pi)`, the whole half-turn). The arrival angles repeat every turn: a bin edge `z` that lies `mq` whole turns beyond the first turn's start brings in `mq` whole turns of photons, each counted with its own `(z - a)`, which is the term `lift` (`mq` the number of whole turns, `Q` the edge brought back into the first turn, with `K(z) = W lift(z) - mq M1 + K(Q)`, `W` the count per turn and `M1` the sum of its arrival angles). Everything is linear in the query's quantities, so the second difference over the bin's two edges and the stretch's two ends is formed first and `G` applied to the differences; all rings, half-turns, nodes and bins make one matrix `Mk`, and the bins are the matrix product `Mw @ Np` with the spectrum at the nodes. A ring that lies whole inside the spot (`dl >= pi`) gives each bin its width times the whole turn's count.

What is exact: the second-difference identity, for the interpolated integrand. What is approximate: the finite number of rings, points and energies (3.4), and the location of `t*`, which is found by inverse interpolation in a table of `a` on an even grid of `ND = 129` points in the Chebyshev variable and not refined. For the exact integrand an error in `t*` changes K only at second order (the summand `z - a` vanishes there); for the interpolated integrand there is also a first-order term, the interpolation's residual at that point times the error. Measured by varying the grid alone from 65 to 2,049 points (`2-python-solution/tools/arrival_grid_check.py`, in `check.py --slow`: "the even grid that locates arrival times: 129 points against 2049"): the working 129 points change the bins by 3.5e-9 of the mean bin and the quick example's ln L at the working setting by 1.8e-6 against the 2,049-point grid; 65 points by 2.4e-8 and 1.3e-5 (the docstring of `fast_kernel.py` gives 4e-9 and 2e-6).

Function: `fast_kernel.spot_counts` (`Q`, `mq`, `lift`, `G`, `q0`, `Mk`); `spot_counts` in `nicer_fast.c`. Checks: those of 3.1 (the plain program has no K and no stretches: it resamples each patch's flux at the arrival phases), and in `check.py` "quick example: ours at the working setting against ours converged" (10 rings per panel, 16 points and 24 energies against 48, 32 and 48; limit 0.02 in ln L).

### 3.4 What is approximate, and what the program refuses

- Rings: Gauss points in latitude, `n_rings` per panel (10 at the working setting, 48 converged). Points: `n` Chebyshev directions per ring and half-turn (16; 32). Energies: `nE` Chebyshev nodes in `ln E` (24; 48; section 4.2). The measured effect of the settings on the quick example's ln L is the table in section 4.3 and the check "quick example: ours at the working setting against ours converged".
- The approximation of the path integrals: 8.5e-13 rad and 1.6e-12 R/c (2.5). Where the integrals themselves are used, 32 Gauss points against 96 change the quick example's ln L by 7e-12 and the second check star's by 1e-6; 3 secant steps for the ring ends against 8 change them by 7e-12 and 2e-6 (the docstring of `fast_kernel.spot_counts`; these two stand in the docstring only, with no committed command).
- The even grid that locates `t*`: 3.5e-9 of the mean bin, 1.8e-6 in ln L (3.3).
- The visible part of each ring's turn is taken as one stretch of directions; where it is in two pieces the program sets the hidden middle to zero at its points, which is approximate, and it has no bound on the error this leaves in a bin. What the committed test finds (13 such rings of 14,968 on 400 random stars over a range wider than the public example's prior, one of them with all its light in the second piece, on a star the compiled program refuses) is in `COVERAGE.md`, section 3.5.
- Covered: `L / R` between 0.12 and 0.62 (`U_RANGE`), and one path of light to every point of the spot. The test: no visible ray sweeps more than the one that leaves edge-on; if that sweep reaches `2 pi` minus the angle of a ring's farthest point, a second path is possible and the program stops with a message (`TwoPaths`) instead of giving a number. For a round star this can first happen at `L / R = 0.568`, and only for spots that reach the far side. The public example's prior stays below 0.4725 and is never refused.
- Whether a ray, once it has left the surface, is caught again by the bulging equator is not tested.

## 4. The likelihood as the public code defines it

This section is the definition of the number to be reproduced. It is not physics of ours: interpolation, sums, one integral per channel, as the public code's source has them. `public_likelihood.py` holds this in our own Python; nothing of the public code (X-PSI) is imported or run.

### 4.1 The definition

Plain: for each energy channel, the chance of the observed counts in the 32 slices of the turn, given the model's counts plus an unknown steady background, averaged over the background between two limits that the public code itself chooses; the logarithm of that, added over the channels.

Exact. For channel r with counts `d_j` (j over the n = 32 phase bins), model source rate `s_j` (the mean count rate in bin j, 4.2), exposure time T and `SCALE = T / n`, the expected counts are `c_j(B) = SCALE (s_j + B)` with B the background count rate, the same in every bin of the channel. The channel's contribution is

    lnL_r = ln INT_{lower}^{upper} PROD_j [ c_j(B)^{d_j} exp(-c_j(B)) / d_j! ] dB,

and the returned number is `sum_r lnL_r` over the 291 channels (the public code's source, spec sections 9.2 to 9.7; `marginal_loglike_channel` and `marginal_loglike` in `public_likelihood.py`; `marginal_loglike` in `nicer_fast.c`). There is no division by the width of the limits. The limits: the public code finds the background rate that makes the product largest, by Newton's steps inside a bracket, starting from the data's rate minus the model's and stopping when the step is below `epsilon = 1e-3` of the width `sigma_B = 1 / sqrt(sum_j d_j / (s_j + B)^2)`; then `lower = B - 10 sigma_B` (not below 0), `upper = B + 10 sigma_B`, both cut to the channel's support (for the example `0.9` to `1.1` times a given background; a negative upper support means no upper limit); if the whole range lies outside the support, the whole support is used (`channel_limits` in `public_likelihood.py` and `nicer_fast.c`). A channel whose model counts exceed the data by more than 20 times their square root ends the whole evaluation (the early exit `slim`). The public code integrates with the GSL routine cquad to a relative tolerance of 1e-8; our reference form uses Gauss-Legendre panels doubled until two answers agree to 1e-12 (`marginal_loglike_channel`; `channel_reference` in C), and our fast form does all channels at once with one Gauss rule on the part of the range where the integrand is above `e^-36` of its largest value, the cuts checked against the true logarithm (`marginal_loglike_fast`). The iterate at which the search stops is kept as the public code keeps it, because the limits are built from it.

### 4.2 The public code's two interpolation steps, and our two routes

The public code works its signal out at 64 energies, spaced evenly in `log10(E)` from the first to the last instrument edge, and at 64 phases from 0 to 1. Energy (spec 8.3): for each phase it forms `y_j = E_j F(E_j) ln(10)` against `x_j = log10(E_j)`, puts the GSL interpolant "akima-periodic" through the 64 points, and takes the exact integral of that interpolant over each of the 291 instrument input intervals (the upper limit cut to the last energy). Phase (spec 8.5): for each channel it puts the same kind of interpolant through the 64 phases, adds phase_shift to each data bin's edges (with a plus sign), brings them into `[0, 1)`, integrates the interpolant over the bin (in two pieces when the bin wraps), adds a piece only if it is positive, sets a negative total to zero, and multiplies by n = 32, so that `s_j` is the mean rate in the bin. The akima-periodic algorithm is written out from GSL's source (`AkimaPeriodic`, `energy_integrate`, `phase_integrate` in `public_likelihood.py`); it wraps the end slopes round whether or not the data are periodic, and the public code uses it for the energy integration too, where they are not.

Our two routes (`quick_example.py`):

- "ours": the bin counts come straight from `fast_kernel.spot_counts` (exact phase-bin integrals, 3.3) at `nE` Chebyshev nodes in `ln E`; the integral of `F(E)` over each instrument interval is `INT (E F) d(ln E)`, taken as the exact integral of the polynomial of degree `nE - 1` in `ln E` through `E F` at the nodes, which is one fixed matrix `W` (`quick_example.cheb_energy_matrix`); then the instrument as the example has it (an area of 1800 cm^2, each channel reading one input interval, as `quick_example.ours` takes it), times n, then the likelihood of 4.1. Neither interpolation step of the public code is applied.
- "emulated": our converged curve sampled at the public code's 64 phases and 64 energies, then the public code's two interpolation steps as written out above, then the likelihood. It shows what those two steps cost in accuracy.

### 4.3 How ours was checked against the public code's own shipped values

Plain: the public code ships four kinds of numbers that depend on no run of it here: two values of the likelihood step alone, the expected likelihood of its quick example with the tolerance its own check accepts, and the expected likelihood of a second check. Ours is compared with each.

- The likelihood step alone: the public code's two shipped unit values, for one channel of 32 counts with a given signal, the background limited to `[0, 1e-11]` and with no upper limit: -4521.763801536406 and -2772.718471175409. Checks in `check.py`: "likelihood step against the public code's shipped unit value ... (relative)", one per value; limit 1e-8, relative. Measured: ours -4521.763801535126 (relative difference 2.8e-13) and -2772.718471175409 (0).
- The quick example (a star of 1.4 Suns and 12 km turning 314 times a second, seen 60 degrees from its axis from 1 kpc, one spot of radius 0.75 rad centred 70 degrees from the axis at 10^6.7 K; 291 channels, 32 bins, 1000 s): the shipped likelihood -32552.616905, which the public code's own check accepts to a relative tolerance of 1e-5, that is 0.3255. Check: "quick example: ours (converged) against the public code's shipped likelihood; its own tolerance" (48 rings per panel, 32 points, 48 energies). Measured: ours converged -32552.710, within the tolerance, 0.093 below the shipped value; at the working setting (10, 16, 24) -32552.719; the emulated route -32552.669, 0.052 below. Which of the two, ours or the shipped one, is nearer the exact value of the model is not shown by this score. The public code, run here at finer settings on 23 November 2026, answers it bin by bin: its counts move toward ours as its mesh and its energies are refined, and ours do not move (`CLAIMS.md`, H8 and H11).
- The second check point (the same data; the star with radius 10 km and 10^6.8 K; no upper limit on the background; the public code at 128 energies): the shipped value -47881.278, tolerance 0.479. Ours converged is -47874.7600, that is +6.518 above the shipped value: ours FAILS this check. Check in `check.py --slow`: "the public code's second check: ours converged minus its shipped value, against the +6.518 found before (ours FAILS that check)" (`second_check_point.py`; it tests that the failure is the same one each time, limit 5e-3 from +6.518).
- The fast likelihood form against the step-by-step form: "quick example: the fast likelihood form against the step-by-step form" (limit 1e-8) and the slow "the fast likelihood form against the step-by-step form on 247 hard made-up channels, per channel" (limit 1e-8; `fast_likelihood_stress.py`: few counts, source far above or below the data, wide and narrow limits). Measured: 2.9e-11 on the quick example and 5.8e-11 per channel on the 247 channels; `fast_likelihood_check.py` prints both check points and three more cases, all within 1.1e-11 in the total.

## 5. Steps not yet written out

1. Section 1, the inputs: that the paths outside the flattened, spinning star are those of P1 to P4 for a round body that does not spin, with the same L, from each point's own radius `R(theta)`. The programs assume it; the page does not derive it.

Everything else on this page is derived from P1 to P4 and H2 with its steps shown, or is an input of the public example as flagged, or is a definition of the public code's own number (section 4). P4 and the arrival of every photon remain IMPORTS wherever they are used.

## 6. What the agreement with the public code does and does not show

What it shows. Our likelihood for the quick example passes the public code's own check by its stated tolerance, at every setting from 8 rings per panel upward, and is steady to about 0.0022 from 12 rings per panel upward. Two programs of ours, which add the light up in different ways (section 3 against small patches with the lensing factor written out), agree to 0.001 in the likelihood and to 4 parts in a million of the mean bin in every bin. They share the path integrals and the star's formulas (shape, tilt, area, energy shift, Doppler factor, the factors of the count formula), so their agreement checks the adding-up and not those shared formulas; the agreed reference curves of the test article and the example's shipped counts file check those.

What it does not show. That every one of the 9,312 bin predictions equals the public code's: a single likelihood cannot show that. Nor, from the score alone, which of the two likelihoods, ours or the shipped one, is nearer the exact value of the model: the shipped one is 0.093 from both of our programs. The public code, run at finer settings on 23 November 2026, answers it bin by bin (`CLAIMS.md`, H8 and H11).

What it means for the physics is said once, in section 1. One more fact belongs beside it: general relativity's premise that the frames share one event horizon is used nowhere in this calculation, and nothing here tests it.
