"""Many stars at once from Python, through the shared library libnicer_fast.so (make builds it), with Python's standard
library only (ctypes), so that a sampler written in Python can ask for many stars in one call. Runs on the Linux side.

    import nicer_fast_lib
    lib = nicer_fast_lib.Batch("tables.bin", 10, 16, 24)          # the tables, then rings per panel, points, energies
    lnl, codes = lib.work_out(stars, threads=8)                   # stars: a list of stars, each eight numbers

The eight numbers of a star and the meaning of the values and codes are those of the program's batch mode (README.md):
-1e300 stands for minus infinity (the public code's early exit, or a refused star, whose code is not 0). The tables and
the settings belong to the process, not to one Batch: the last Batch made, or the last choose, sets them for all.

usage as a command: python3 -I -B nicer_fast_lib.py <tables> <rings per panel> <points> <energies> <threads> <stars file>
prints one line per star exactly as the program's batch mode prints it (check.py compares the two)."""
import ctypes
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


class Batch:
    def __init__(self, tables, npp, n, nE, library=os.path.join(HERE, "libnicer_fast.so")):
        self.lib = ctypes.CDLL(library)
        self.lib.nicer_load.argtypes = [ctypes.c_char_p]
        self.lib.nicer_choose.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        self.lib.nicer_choose.restype = None
        self.lib.nicer_batch.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)]
        self.lib.nicer_load(os.fsencode(tables))       # once per process: a second Batch keeps the tables of the first
        self.lib.nicer_choose(npp, n, nE)

    def choose(self, npp, n, nE):
        self.lib.nicer_choose(npp, n, nE)

    def work_out(self, stars, threads=1):
        count = len(stars)
        p = (ctypes.c_double * (8 * count))(*[float(x) for star in stars for x in star])
        lnl = (ctypes.c_double * count)()
        codes = (ctypes.c_int * count)()
        if self.lib.nicer_batch(count, p, threads, lnl, codes, None):
            raise MemoryError("the library could not have a thread's working memory")
        return list(lnl), list(codes)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(__doc__.strip().splitlines()[-2])
        sys.exit(2)
    tables, npp, n, nE, threads, path = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    with open(path, encoding="ascii") as fh:
        stars = [[float(x) for x in line.split()] for line in fh if line.strip()]
    lnl, codes = Batch(tables, npp, n, nE).work_out(stars, threads)
    for v, c in zip(lnl, codes):
        print("-inf %d" % c if v < -1.0e299 else "%s %d" % (format(v, ".17g"), c))
