# -*- coding: utf-8 -*-
"""Writes the tables that the compiled program (c/nicer_fast.c) reads: everything that is made once and never changes
between workings-out. All of it comes from our own Python programs and from the public example's data files:
  the two-variable fit of light's path (fast_kernel._path_fit), the table for the latitude cuts (fast_kernel._tables),
  Gauss-Legendre points, the Chebyshev tools for the chosen numbers of points (fast_kernel._cheb),
  the energy points and the matrix of energy integrals (quick_example.cheb_energy_matrix),
  the public example's counts, its background limits and the factorial sums, and the public code's constants.
File format: for each array, int32 length of the name, the name, int32 number of dimensions, int32 sizes, float64 values.

usage: python export_for_c.py <the example's Data folder> <output file>
"""
import math
import os
import struct
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
import numpy as np
import fast_kernel as fk
import public_likelihood as pl
import quick_example as qe
from consts import BB_CONST, C_CM_S, GM_SUN

data_dir, out_path = sys.argv[1], sys.argv[2]
arrays = []


def put(name, a):
    arrays.append((name, np.ascontiguousarray(np.asarray(a, dtype=np.float64))))


Ct, Ca, ZM, ku, kz = fk._path_fit()
put("fit_ct", Ct)
put("fit_ca", Ca)
put("fit_meta", [fk.FIT_U[0], fk.FIT_U[1], ZM, fk.FIT_XMIN, fk.FIT_TILT])
ug, t0, t0x, _, t0xx = fk._tables()
put("vis_ug", ug)
put("vis_t0", t0)
put("vis_t0x", t0x)
put("vis_t0xx", t0xx)
for npts in (6, 8, 10, 12, 14, 16, 20, 24, 28, 40, 48, 64):
    g, w = np.polynomial.legendre.leggauss(npts)
    put("gl_g_%d" % npts, 0.5 * (g + 1.0))
    put("gl_w_%d" % npts, 0.5 * w)
for n in (12, 16, 24, 32):
    tj, costj, sintj, G, q0, xid, to_d, kk = fk._cheb(n)
    put("ch_cos_%d" % n, costj)
    put("ch_sin_%d" % n, sintj)
    put("ch_G_%d" % n, G)                # (n + 1, n)
    put("ch_q0_%d" % n, q0)
    put("ch_tod_%d" % n, to_d)           # (n, ND)
put("ch_meta", [fk.ND])
edges_E = pl.example_instrument_edges()
for nE in (12, 16, 24, 32, 48):
    E, W = qe.cheb_energy_matrix(edges_E[0], edges_E[-1], nE, edges_E)
    put("en_E_%d" % nE, E)
    put("en_W_%d" % nE, qe.AREA * qe.NB * W)         # (channels, energies): mean count rate in a bin per unit photon flux
counts, bg, support = qe.load_data(data_dir)
put("data_counts", counts)
put("data_support", support)
put("data_fact", pl.log_factorial_sums(counts))
# constants: light speed (cm/s), G M_sun / c^2 (cm) as the public code has it, kpc (cm), k_B / keV, the black-body constant,
# exposure (s), spin (Hz), number of phase bins
put("consts", [C_CM_S, qe.X_GM_M * 100.0, qe.X_KPC_M * 100.0, qe.X_KB_OVER_KEV, qe.X_BB, qe.T_EXP, qe.P["freq"], qe.NB])
with open(out_path, "wb") as fh:
    fh.write(struct.pack("<i", len(arrays)))
    for name, a in arrays:
        nb = name.encode("ascii")
        fh.write(struct.pack("<i", len(nb)))
        fh.write(nb)
        fh.write(struct.pack("<i", a.ndim))
        for s in a.shape:
            fh.write(struct.pack("<i", s))
        fh.write(a.tobytes(order="C"))
print("wrote %d arrays, %d bytes: %s" % (len(arrays), os.path.getsize(out_path), out_path))
