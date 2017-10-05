from fpylll import IntegerMatrix, LLL
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.bkz_stats import dummy_tracer
from fpylll import Enumeration, EnumerationError 
from fpylll.util import gaussian_heuristic


def load_challenge_and_randomize(n):
    A_pre = IntegerMatrix.from_file("svpchallenge/svpchallengedim%dseed0.txt"%n)
    LLL.reduction(A_pre)
    A = IntegerMatrix.from_matrix(A_pre, int_type="long")
    params = fplll_bkz.Param(block_size=n, max_loops=1, strategies=fplll_bkz.DEFAULT_STRATEGY, flags=fplll_bkz.GH_BND, min_success_probability=.8)
    bkz = BKZReduction(A)
    bkz.lll_obj()
    bkz.randomize_block(0, n, density=n/4)
    bkz.lll_obj()   
    return A, bkz

## Re-implement bkz2.svp_reduction, with a precise radius goal rather than success proba

def svp_reduction_until_goal(bkz, params, goal):
    n = bkz.M.d
    r = [bkz.M.get_r(i, i) for i in range(0, n)]
    gh = gaussian_heuristic(r)

    while bkz.M.get_r(0, 0) > goal:
        bkz.randomize_block(0, n)
        bkz.svp_preprocessing(0, n, params)

        strategy = params.strategies[n]
        radius = goal
        pruning = strategy.get_pruning(goal, gh)

        try:
            enum_obj = Enumeration(bkz.M)
            max_dist, solution = enum_obj.enumerate(0, n, radius, 0, pruning=pruning.coefficients)[0]
            bkz.svp_postprocessing(0, n, solution, tracer=dummy_tracer)
            rerandomize = False

        except EnumerationError:
            rerandomize = True

    bkz.lll_obj()
    return 