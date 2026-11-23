"""Stars drawn at random from the public quick example's prior, one star per line, for the batch mode of nicer_fast.

The ranges are those of 2-python-solution/tools/prior_draws_c_and_python.py: mass 1.0 to 1.6, radius 10 to 13, distance 0.5
to 2, cos i 0 to 1, phase shift -0.25 to 0.75, spot colatitude and radius 0.001 to pi/2 - 0.001, log10 T 6 to 7, drawn in that
order. The random numbers are Python's own (the standard library's random module, whose generator gives the same numbers on
every machine for the same seed), so this file runs on the Linux side, whose Python may have no NumPy; the sixty stars of
prior_draws_c_and_python.py come from NumPy's generator and are not these.
usage: python3 -I -B draw_stars.py <seed> <number of stars> <file>
Each number is written as Python's repr, which the program reads back to the same bits."""
import math
import random
import sys

if len(sys.argv) != 4:
    print(__doc__.strip().splitlines()[-2])
    sys.exit(2)
seed, count, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = random.Random(seed)
with open(path, "w", encoding="ascii", newline="\n") as fh:
    for k in range(count):
        p = [rng.uniform(1.0, 1.6), rng.uniform(10.0, 13.0), rng.uniform(0.5, 2.0), rng.uniform(0.0, 1.0), rng.uniform(-0.25, 0.75),
             rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(0.001, math.pi / 2 - 0.001), rng.uniform(6.0, 7.0)]
        fh.write(" ".join(repr(x) for x in p) + "\n")
print("%d stars from the public example's prior (seed %d): %s" % (count, seed, path))
