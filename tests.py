import itertools
import math
import numpy as np
import scipy.stats as stats
import unittest

import networkx as nx
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

def assert_compnents(g, n=1):
    ccs = nx.connected_components(g)
    n_ccs = len(ccs)
    if n != n_ccs:
        raise AssertionError("Number of connected components wrong: %s!=%s (%s)"%(
            n, n_ccs, [len(cc) for cc in ccs]))


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


class _TestRandom(unittest.TestCase):
    model_name = None

    def test_size(self):
        M = bm.get_model(self.model_name)
        assert len(M.t(0)) == 128
    def test_all_managed(self):
        assert_all_managed(bm.get_model(self.model_name))
        assert_all_managed(bm.get_model(self.model_name, q=8))
        assert_all_managed(bm.get_model(self.model_name, p_in=.5, p_out=.5))


    def test_random(self):
        #class Binom():
        #    def __init__(self, n, p):
        #        self.n, self.p = n, p
        #    def rvs(self):
        #        return int(round(self.n*self.p))
        #import scipy.stats
        #scipy.stats.binom = Binom

        def getM():
            return bm.get_model(self.model_name, p_in=.5, p_out=.5)
        assert_distribution(lambda: getM().t(0).number_of_edges(),
                            128*127/2 * .5, std=None, p=.01)

        def getM():
            return bm.get_model(self.model_name, p_in=.6, p_out=.6)
        assert_distribution(lambda: getM().t(0).number_of_edges(),
                            128*127/2 * .6, std=None, p=.01)

        def getM():
            return bm.get_model(self.model_name, p_in=.6, p_out=.6, q=8)
        assert_distribution(lambda: getM().t(0).number_of_edges(),
                            256*255/2 * .6, std=None, p=.01)

        def getM():
            return bm.get_model(self.model_name, p_in=.6, p_out=.6)
        assert_distribution(lambda: getM().t(50).number_of_edges(),
                            128*127/2 * .6, std=None, p=.01)

        def getM():
            return bm.get_model(self.model_name, p_in=.6, p_out=.6)
        assert_distribution(lambda: getM().t(75).number_of_edges(),
                            128*127/2 * .6, std=None, p=.01)

        def getM():
            return bm.get_model(self.model_name, p_in=.6, p_out=.6)
        assert_distribution(lambda: getM().t(25).number_of_edges(),
                            128*127/2 * .6, std=None, p=.01)


    def test_ccs(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        g = M.t(0)
        assert_compnents(g, getattr(self, 'n_ccs', 4))

        g = M.t(50)
        assert_compnents(g, getattr(self, 'n_ccs', 4))

        g = M.t(25)
        assert_compnents(g, getattr(self, 'n_ccs_off', 4))

    def test_num_comms(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        assert_equal(len(M.comms(0)),   4)
        assert_equal(len(M.comms(25)),  4)
        assert_equal(len(M.comms(50)),  4)
        assert_equal(len(M.comms(85)),  4)
        assert_equal(len(M.comms(100)), 4)

class TestMerge(_TestRandom):
    model_name = 'StdMerge'
    n_ccs = 3  # 3 conected components
    n_ccs_off = 2
    def test_ccs(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        g = M.t(0);   assert_compnents(g, 3)
        g = M.t(50);  assert_compnents(g, 3)
        g = M.t(25);  assert_compnents(g, 2)
    def test_num_comms(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        assert_equal(len(M.comms(0)),   3)
        assert_equal(len(M.comms(25)),  4)
        assert_equal(len(M.comms(50)),  3)
        assert_equal(len(M.comms(85)),  4)
        assert_equal(len(M.comms(100)), 3)
class TestGrow(_TestRandom):
    model_name = 'StdGrow'
    def test_ccs(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        g = M.t(0);   assert_compnents(g, 4)
        g = M.t(50);  assert_compnents(g, 4)
        g = M.t(25);  assert_compnents(g, 4)
class TestMixed(_TestRandom):
    model_name = 'StdMixed'
    def test_ccs(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        g = M.t(0);   assert_compnents(g, 4)
        g = M.t(50);  assert_compnents(g, 3)
        g = M.t(25);  assert_compnents(g, 3)
    def test_num_comms(self):
        M = bm.get_model(self.model_name, p_in=.5, p_out=0)
        assert_equal(len(M.comms(0)),   4)
        assert_equal(len(M.comms(25)),  4)
        assert_equal(len(M.comms(50)),  3)
        assert_equal(len(M.comms(85)),  4)
        assert_equal(len(M.comms(100)), 4)

