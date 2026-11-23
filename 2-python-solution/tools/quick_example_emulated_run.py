import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import quick_example as qe
import public_likelihood as pl
counts, bg, support = qe.load_data(sys.argv[1])
edges_E = pl.example_instrument_edges()
for nr, n, nd in ((64, 24, 769), (64, 32, 3073)):
    t0 = time.perf_counter()
    v = qe.emulated(qe.P, counts, support, edges_E, n_rings=nr, n=n, nd=nd)
    print("the public code's two interpolation steps applied to our curve (%d rings per panel, %d points, grid %d): ln L = %.6f ; minus expected: %+.6f  (%.0f s)"
          % (nr, n, nd, v, v - qe.EXPECTED, time.perf_counter() - t0), flush=True)
