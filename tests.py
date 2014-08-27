import itertools
import math
import numpy as np
import scipy.stats as stats

from nose.tools import *


def assert_distribution(f, mean, std, p=.1):
    vals = [ ]

    for N in [10, 20, 50, 100, 500]:
        while len(vals) < N:
            vals.append(f())
        #print vals

        smean, sstd = np.mean(vals), np.std(vals)
        if sstd == 0:
            sstd = 1e-6

        t = (mean - smean) / float(sstd / math.sqrt(len(vals)))
        t = -abs(t)
        print t, mean, smean, sstd

        p1 = stats.t.cdf(t, len(vals)-1)
        assert p1 == abs(p1)

        p2 = 2*p1
        print "p:", p2

        if p2 > p:
            print "passed with N=%s"%N
            return True   # test passed
    raise AssertionError("Test failed, %s!=%s(std=%s) (%s !< %s) with %s samples"%(mean, smean, sstd, p, p2, N))

def test_assert_dist():
    assert_distribution(stats.norm(0).rvs, 0, 1)


def test():
    f = lambda : 1

    assert_distribution(f, 1, 0)

import bm

def assert_all_managed(bm):
    """Check each pair of nodes and make sure that its link is managed.

    Iterate through all pairs of nodes, and for all managers of the
    benchmark object, make sure that exactly one manager claims to be
    managing its link."""
    for a, b in itertools.combinations(bm.g.nodes_iter(), 2):
        found = False
        managed_by = [mgr.manages(a, b) for mgr in bm.managers]
        # True sums to one, False to zero.  Each node should be
        # managed by exactly one manager.
        if sum(managed_by) != 1:
            raise AssertionError('Nodes %s and %s are not managed (%s)'%(
                a, b, managed_by))

def test_all_managed():
    assert_all_managed(bm.get_model('StdGrow' ))
    assert_all_managed(bm.get_model('StdMerge'))
    assert_all_managed(bm.get_model('StdMixed'))

    assert_all_managed(bm.get_model('StdGrow' , q=8))
    assert_all_managed(bm.get_model('StdMerge', q=8))
    assert_all_managed(bm.get_model('StdMixed', q=8))


def test_1():
    M = bm.get_model('StdGrow')
    assert len(M.t(0)) == 128
