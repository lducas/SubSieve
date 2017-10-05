from numpy import zeros, float64, int64, matrix, array
from fpylll.util import gaussian_heuristic
from time import time
import ctypes
import _ctypes
from math import floor, log
import subprocess
import threading


def c_double_ptr(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def c_long_ptr(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


class SubSieveLib(object):
    def __init__(self, n, d, gso):
        self.n = n
        self.lib = ctypes.cdll.LoadLibrary("./SubSieveLib.so")
        self.d = d
        self.r = [gso.get_r(i, i) for i in range(self.n)]
        self.gh = gaussian_heuristic(self.r[d:])
        self.gs = zeros((self.n, self.n), dtype=float64)
        for i in xrange(self.n):
            self.gs[i][i] = gso.get_r(i, i)
            for j in xrange(i):
                self.gs[i][j] = gso.get_mu(i, j)

        self.lib.initialize(self.n, self.d, c_double_ptr(self.gs), ctypes.c_double(self.gh))

    def sieve(self):
        t = threading.Thread(target=self.lib.sieve)
        t.daemon = True
        t.start()
        while t.is_alive(): # wait for the thread to exit
            t.join(.01)

    def lift(self):
        self.gh_full = gaussian_heuristic(self.r)
        self.lib.lift(self.n, self.d, c_double_ptr(self.gs), ctypes.c_double(self.gh_full))

    def insert_projector(self, v):
        vv = zeros(self.n, dtype=int64)
        for i in range(self.n):
            vv[i] = v[i]
        self.lib.insert_projector(c_long_ptr(vv))

    def export_db(self, nb):
        db_e = zeros((nb, self.n), dtype=int64)
        self.lib.export_db(nb, c_double_ptr(db_e))
        return db_e

    def import_db(self, db_e):
        (nb, _) = db_e.shape
        self.lib.import_db(nb, c_double_ptr(db_e))

    def db_size(self):
        return self.lib.db_size()

    def resize_db(self, N):
        self.lib.resize_db(N)

    def __del__(self):
        _ctypes.dlclose(self.lib._handle)
