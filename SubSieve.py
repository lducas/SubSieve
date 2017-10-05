# -*- coding: utf-8 -*-

from fpylll import IntegerMatrix, LLL, GSO
from fpylll import BKZ as fplll_bkz
from fpylll.util import gaussian_heuristic
from middleware import SubSieveLib
from random import randint
from numpy import int64


def update_basis(lll, V):
    M = lll.M
    l = len(V)
    n = M.d

    for v in V:
        M.create_row()
        with M.row_ops(M.d-1, M.d):
            for i in range(n):
                M.row_addmul(M.d-1, i, v[i])    

    for k in range(l):
        M.move_row(M.d-1, 0)

    lll()

    for i in range(l):
        M.remove_last_row()


def extract_good_subbasis(siever, m):
    B = []
    N0 = siever.db_size()
    N = N0
    for i in range(m):
        v = siever.export_db(1)[0]
        B += [v]        
        N = max(3*N/4, N0/m)        
        siever.resize_db(N)
        siever.insert_projector(v)

    return B


def sub_sieve_plus(lll, d):

    M = lll.M
    n = M.d

    siever = SubSieveLib(n, d, M)
    siever.sieve()
    siever.lift()

    subB = extract_good_subbasis(siever, d)
    update_basis(lll, subB)
    del siever


def iterated_sub_sieve(A, goal):
    n = A.nrows
    M = GSO.Mat(A)
    lll = LLL.Reduction(M)
    lll()
    n = M.d
    d = n/4
    while M.get_r(0, 0) > goal:
        sub_sieve_plus(lll, d)
        d -= 1
    return d+1