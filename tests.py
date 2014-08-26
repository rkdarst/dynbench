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
    raise AssertionError("Test failed, (%s !> %s) with %s samples"%(p2, p, N))

def test_assert_dist():
    assert_distribution(stats.norm(0).rvs, 0, 1)


def test():
    f = lambda : 1

    assert_distribution(f, 1, 0)
