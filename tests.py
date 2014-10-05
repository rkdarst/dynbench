import itertools
import math
import numpy as np
import os
import scipy.stats as stats
import shutil
import tempfile
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


    def test_output_bynode(self):
        """Ensure that main writes output: bynode"""
        tmpdir = tempfile.mkdtemp(prefix='dynbench-')
        try:
            M = bm.main(argv=[bm.__file__, self.model_name,
                              tmpdir+'/output-1', # outdir
                              #'--comm-format=bynode', # this should be default
                              ])
            assert_equal(len(os.listdir(tmpdir)), 202, "202 files should be written")
            comms = M.comms(2)
            for line in open(tmpdir+'/output-1.t00002.comms'):
                if line.startswith("#"): continue
                n, c = line.split()
                assert int(n) in comms[int(c)]
        finally:
            shutil.rmtree(tmpdir)

    def test_output_oneline(self):
        """Ensure that main writes output: oneline"""
        tmpdir = tempfile.mkdtemp(prefix='dynbench-')
        try:
            M = bm.main(argv=[bm.__file__, self.model_name,
                              tmpdir+'/output-1', # outdir
                              '--comm-format=oneline',
                              ])
            assert_equal(len(os.listdir(tmpdir)), 202, "202 files should be written")
            comms = M.comms(2)
            for line in open(tmpdir+'/output-1.t00002.comms'):
                if line.startswith('# label: '):
                    c = line[9:]
                if line.startswith("#"): continue
                nodes = line.split()
                for n in nodes:
                    assert int(n) in comms[int(c)]
        finally:
            shutil.rmtree(tmpdir)

    def test_seed(self):
        """Test that seeding does """
        argv = [bm.__file__, self.model_name,
                '--p_in=.5', '--p_out=.1']
        import random
        random.seed(10)
        M1a, args = bm.main_argv(argv=argv+['--seed=51'])
        random.seed(10)
        M2,  args = bm.main_argv(argv=argv+['--seed=965'])
        random.seed(10)
        M1b, args = bm.main_argv(argv=argv+['--seed=51'])

        def edgeset(g):
            print len(g), g.number_of_edges()
            return set(frozenset((a, b)) for a,b in g.edges_iter())

        assert_equal(    edgeset(M1a.t(5)), edgeset(M1b.t(5)))
        #assert_not_equal(edgeset(M2.t(5) ), edgeset(M1a.t(5)))
        #assert_not_equal(edgeset(M2.t(5) ), edgeset(M1b.t(5)))




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
        assert_equal(len(M.comms(75)),  4)
        assert_equal(len(M.comms(100)), 3)
    def test_det_limit(self):
        Mnd = bm.get_model(self.model_name, p_in=.5, p_out=0,
                           opts=dict(no_det_limit=True))
        Md = bm.get_model(self.model_name, p_in=.5, p_out=0,
                          opts=dict(no_det_limit=False))
        # Test the command line no-det-limit option works.
        argv = [bm.__file__, self.model_name, '/nonexistant', # outdir
                '--graph-format=null','--comm-format=null']
        Mnd2, args = bm.main_argv(argv=argv+['--no-det-limit'])
        # Default shourd be detectability limit.
        Md2, args = bm.main_argv(argv=argv)
        assert_equal(len(Mnd.comms(0)),  3)
        assert_equal(len(Mnd.comms(1)),  4)

        assert_equal(len(Md.comms(0)),   3)
        assert_equal(len(Md.comms(1)),   3)

        assert_equal(len(Mnd2.comms(0)),  3)
        assert_equal(len(Mnd2.comms(1)),  4)

        assert_equal(len(Md2.comms(0)),  3)
        assert_equal(len(Md2.comms(1)),  3)

        #tmpdir = tempfile.mkdtemp(prefix='dynbench-')
        #try:


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

