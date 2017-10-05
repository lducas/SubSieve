# -*- coding: utf-8 -*-

from fpylll import IntegerMatrix, LLL
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.util import gaussian_heuristic

from SubSieve import iterated_sub_sieve
from util import load_challenge_and_randomize, svp_reduction_until_goal

from time import time

prelim_rep = 10
sample_rep = 20

## Re-implement bkz2.svp_reduction, with a precise radius goal rather than success proba

START = time()
for n in range(60, 130, 2):

    AV_ENUM = 0
    AV_SIEVE = 0
    params = fplll_bkz.Param(block_size=n, max_loops=1, strategies=fplll_bkz.DEFAULT_STRATEGY, flags=fplll_bkz.GH_BND)

    print "---------------------", n

    A, bkz = load_challenge_and_randomize(n)
    gh = gaussian_heuristic([bkz.M.get_r(i, i) for i in range(n)])
    goal = None

    for x in range(prelim_rep):
        A, bkz = load_challenge_and_randomize(n)
        
        ENUM_START = time()
        bkz.svp_reduction(0, n, params)
        ENUM_TIME = time() - ENUM_START

        r0 = bkz.M.get_r(0, 0)
        if goal is None:
            goal = 1.002 * r0
        else:
            goal = min(goal, 1.001 * r0)
        print "PREL_ENUM  ", {"n": n, "Time": round(ENUM_TIME, 2), "ghf": round(r0/gh, 4)}

    print "goal length set to", round(goal/gh, 4), "* gh"
    print

    for x in range(sample_rep):
        ########### Timing enum
 
        A, bkz = load_challenge_and_randomize(n)

        ENUM_START = time()
        svp_reduction_until_goal(bkz, params, goal)
        ENUM_TIME = time() - ENUM_START

        r0 = bkz.M.get_r(0, 0)
        print "OPPORT_ENUM", {"n": n, "Time": round(ENUM_TIME, 2), "ghf": round(r0/gh, 4)}
        AV_ENUM += ENUM_TIME / sample_rep

        ########### Timing SubSieve
        A, bkz = load_challenge_and_randomize(n)
       
        SUBSIEVE_START = time()
        d = iterated_sub_sieve(A, goal=goal)
        SUBSIEVE_TIME = time() - SUBSIEVE_START
        
        bkz = BKZReduction(A)
        bkz.lll_obj()
        r0 = bkz.M.get_r(0, 0)

        print "SUB_SIEVE  ", {"n": n, "Time": round(SUBSIEVE_TIME, 2), "ghf": round(r0/gh, 4)}, {"d": d}
        AV_SIEVE += SUBSIEVE_TIME / sample_rep

    print
    print "OPPORT_ENUM average :", AV_ENUM
    print "SUB_SIEVE   average :", AV_SIEVE
    print "ratio :", AV_SIEVE / AV_ENUM
