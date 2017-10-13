# -*- coding: utf-8 -*-

from fpylll import IntegerMatrix, LLL, GSO
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.bkz_stats import dummy_tracer
from time import time
from fpylll import Enumeration, EnumerationError
import sys
from fpylll.util import gaussian_heuristic
from middleware import SubSieveLib


for n in range(40, 90, 2):

    A_pre = IntegerMatrix.from_file("svpchallenge/svpchallengedim%dseed0.txt"%n)
    LLL.reduction(A_pre)
    A = IntegerMatrix.from_matrix(A_pre, int_type="long")
    bkz = BKZReduction(A)
    bkz.lll_obj()
    bkz.randomize_block(0, n, density=n/4)
    bkz.lll_obj()
    
    START = time()
    siever = SubSieveLib(n, 0, bkz.lll_obj.M)
    siever.sieve()
    TIME = time() - START
    print (n, TIME)
    del siever

