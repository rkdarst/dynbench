import math
import numpy as np
import scipy.stats as stats

from nose.tools import *


def assert_distribution(f, mean, std):
    vals = [ ]

    for x in range(10):
        vals.append(f())

    smean, sstd = np.mean(vals), np.std(mean)
    if sstd == 0:
        sstd = 1e9

    t = (x - smean) / float(sstd / math.sqrt(len(vals)))
    print t

    p = stats.t.cdf(t, 9)
    assert p == abs(p)

    assert p > .5


def test():
    f = lambda : 1

    assert_distribution(f, 1, 0)
