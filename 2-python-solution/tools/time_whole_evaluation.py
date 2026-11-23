"""Time of one whole likelihood evaluation of the public code's quick example with our programs, on this machine, one
process, one thread: the star's bin counts (fast_kernel), the energy integrals, the likelihood (fast form).
usage: python -I -B -X utf8 time_whole_evaluation.py [profile]        (the public example's data folder comes from settings.cfg;
       run it when nothing else is running on the machine; "profile" adds a profile of the working setting)"""
import math, os, sys, time, cProfile, pstats
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_v] = "1"          # one thread, set before NumPy loads, so that the time is for one processor core
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quick_example as qe
import public_likelihood as pl
import fast_kernel as fk
import settings
from consts import BB_CONST
counts, bg, support = qe.load_data(settings.data_dir())
fact = pl.log_factorial_sums(counts)
edges_E = pl.example_instrument_edges()
ph_edges = np.linspace(0.0, 1.0, qe.NB + 1)
ref = qe.ours(qe.P, counts, support, edges_E, 48, 32, 48)

def make(nr, n, nE):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    Wn = (qe.AREA * qe.NB * (qe.X_BB / BB_CONST)) * W
    def like(p):
        mass, radius, distance, cos_i, shift, colat, rho, logT = p
        M_eff = mass * qe.X_GM_M * 100.0 * qe.C_CM_S**2 / qe.GM_SUN
        S = fk.spot_counts(M_eff, radius, 314.0, distance * qe.X_KPC_M * 100.0, math.acos(cos_i), colat, rho,
                           qe.X_KB_OVER_KEV * 10.0**logT, E, ph_edges, phase0=0.5 - shift, n_rings=nr, n=n)
        return pl.marginal_loglike_fast(counts, Wn @ S.T, qe.T_EXP, support, fact)
    return like

p = [1.4, 12.0, 1.0, math.cos(60 * math.pi / 180), 0.0, 70 * math.pi / 180, 0.75, 6.7]
print("converged value of ours %.6f ; the public code's expected value %.6f (its tolerance 0.3255)" % (ref, qe.EXPECTED))
for nr, n, nE in ((12, 16, 24), (10, 16, 24), (10, 16, 16), (8, 16, 16), (8, 12, 16)):
    like = make(nr, n, nE)
    v = like(p)
    best = 1e9
    for rep in range(5):
        t0 = time.perf_counter()
        for _ in range(40):
            like(p)
        best = min(best, (time.perf_counter() - t0) / 40)
    print("%2d rings per panel, %2d points, %2d energies: ln L = %.6f ; minus ours converged %+.4f ; minus the public code's expected %+.4f ; %.2f ms per evaluation"
          % (nr, n, nE, v, v - ref, v - qe.EXPECTED, best * 1e3), flush=True)
if len(sys.argv) > 1 and sys.argv[1] == "profile":
    like = make(10, 16, 24)
    pr = cProfile.Profile()
    pr.enable()
    for _ in range(200):
        like(p)
    pr.disable()
    st = pstats.Stats(pr)
    st.sort_stats("tottime").print_stats(18)
