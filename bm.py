import contextlib
import math
import random
import re
import sys
import threading
import time

import networkx as nx
import numpy.random
import scipy.stats

import logging
logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

class DisconnectedError(Exception):
    pass


def k_out_limit(kin, q): return kin-.5*( math.sqrt((q-1)**2+(4*q*kin)) - (q-1))


_numpy_rng_lock = threading.Lock()
@contextlib.contextmanager
def override_numpy_seed(rng):
    """Temporarily override the numpy random number generator seed.

    scipy.stats uses numpy.random without an interface to set the seed or random number generator state machine.  This is a problem, since we need seeding support.  One could just set numpy.random"""
    with _numpy_rng_lock:
        # Ggt old state
        old_state = numpy.random.get_state()
        try:
            numpy.random.seed(rng.randint(0, 2**31-1))
            yield
        except:
            numpy.random.set_state(old_state)
            raise
        finally:
            # Restore old state no matter what.
            numpy.random.set_state(old_state)


class Benchmark(object):
    def __init__(self, p_in=1, p_out=0, tau=100, opts={}):
        self.rng = random.Random(opts.get('seed', None))
        self.opts = opts

        n = 32

        self.c1 = c1 = set(range(0, n  ))
        self.c2 = c2 = set(range(n, n*2))

        nodes = c1 | c2

        managers = [#Static(self, c1, p=p_in),
                    #Static(self, c2, p=p_in),
                    #Merging(self, c1,c2, p_low=p_out, p_high=p_in, tau=tau),
                    ExpandContract(self, c1,c2,
                                   p_in=p_in, p_out=p_out, tau=tau),
                    ]
        self.managers = managers
        self.g = g = nx.Graph()

        for c in (c1, c2):
            for n in c:
                g.add_node(n)

    def graph(self, t):
        g = self.g.copy()
        for mgr in self.managers:
            mgr.t(g, t)
        return g
    t = graph
    def comms(self, t):
        comms = { }
        for mrg in self.managers:
            for cname, cnodes in mrg.comms(t).iteritems():
                if cname in comms:
                    raise ValueError("Duplicate community name: %s"%cname)
                comms[cname] = cnodes
        return comms


def shuffled(rng, x):
    if isinstance(x, list):
        x = x[:]    # force a copy
    else:
        x = list(x)
    rng.shuffle(x)
    return x

def add_edge_nonexists(g, n1, n2):
    assert not g.has_edge(n1, n2), "Graph has %s-%s."%(n1, n2)
    g.add_edge(n1, n2)

def choose_random_edges(c1, c2=None, m=None, rng=None):
    # One community internal edges
    if c2 is None:
        one_cmty = True
        n_links = n_links = len(c1) * (len(c1)-1) / 2
        c2 = c1
    # Two communities external edges
    else:
        one_cmty = False
        n_links = len(c1) * len(c2)
        assert len(c1 & c2) == 0, "c1 and c2 overlap, NotImplemented"

    # Sparse links
    if m < .5 * n_links:
        # Choose edges for sparse graphs.
        edges = set()
        for _ in range(m):
            while True:
                n1 = rng.sample(c1, 1)[0]#sets support sample but not choice
                n2 = rng.sample(c2, 1)[0]
                #print n1, n2
                if n1 == n2: continue
                e = frozenset((n1, n2))
                if e in edges: continue
                edges.add(e)
                break
        assert len(edges) == m
        edges = list(edges)
    # Dense links
    else:
        if one_cmty:
            lst = list(c1)
            edges = set(frozenset((lst[i], lst[j]))
                      for i in range(len(lst)) for j in range(i+1, len(lst)) )
            assert len(edges) == len(c1) * (len(c1)-1) / 2
        else:
            edges = set(frozenset((a, b))  for a in c1  for b in c2)
            assert len(edges) == len(c1)*len(c2), "overlap between c1 and c2?"
        possible_edges = list(edges)
        rng.shuffle(possible_edges)

        edges = possible_edges[:m]
        assert len(set(edges)) == m
    return edges



class Static(object):
    def __init__(self, bm, c1, c2=None, p=None):
        self.c1 = c1
        self.c2 = c2
        self.bm = bm
        assert p is not None

        debug("Static, c2=%s", type(c2))
        # Number of total possible edges
        if c2 is None:
            n_links = len(c1) * (len(c1)-1) / 2
        else:
            n_links = len(c1) * len(c2)
        debug("Static, meanlinks=%s, n_links=%s, p=%s", p*n_links, n_links, p)
        # Number of actual edges
        if not self.bm.opts.get('Gnm', False):
            # Gnp random graph ensemble
            with override_numpy_seed(self.bm.rng):
                n_edges = scipy.stats.binom(n_links, p).rvs()
        else:
            # Gnm random graph ensemble
            n_edges = int(round(n_links * p ))

        self.edges_active = choose_random_edges(c1=c1, c2=c2, m=n_edges,
                                           rng=self.bm.rng)
        debug("Static, links=%s", len(self.edges_active))

    def t(self, g, t):
        for a,b in self.edges_active:
            add_edge_nonexists(g, a, b)
        #g.add_edges_from(self.edges_active)

        if self.c2 is None:
            if len(nx.connected_components(g.subgraph(self.c1))) != 1:
                raise DisconnectedError("Subgraph is disconnected (Static object)")
    def comms(self, t):
        return { }
    def manages(self, a, b):
        """Return true if two nodes link is managed by this object"""
        if self.c2 is None:
            if a in self.c1 and b in self.c1:
                return True
        else:
            if (   (a in self.c1 and b in self.c2)
                or (b in self.c1 and a in self.c2)):
                return True
        return False


class Merging(object):
    def __init__(self, bm, c1, c2, p_low, p_high, tau, phasefactor=0.,
                 c_id_1=0, c_id_2=1, c_id_merged=2):
        self.bm = bm
        self.n_links = n_links = len(c1) * len(c2)
        debug("Merging, meanlinks_low=%s, meanlinks_high=%s", p_low*n_links, p_high*n_links)
        self.c1 = c1
        self.c2 = c2
        self.c_id_1 = c_id_1
        self.c_id_2 = c_id_2
        self.c_id_merged = c_id_merged

        self.p_low = p_low
        self.p_high = p_high
        if not self.bm.opts.get('Gnm', False):
            # Gnp random graph ensemble
            with override_numpy_seed(self.bm.rng):
                self.m_low  = scipy.stats.binom(n_links, p_low ).rvs()
                self.m_high = scipy.stats.binom(n_links, p_high).rvs()
        else:
            # Gnm random graph ensemble
            self.m_low  = int(round(n_links * p_low))
            self.m_high = int(round(n_links * p_high))

        debug("Merging, links_low=%s, links_high=%s, p_low=%s, p_high=%s",
              self.m_low, self.m_high, p_low, p_high)
        self.tau = tau
        self.phasefactor = phasefactor
        self.p_limit =  k_out_limit(self.p_high*len(c1), 2) / float(len(c1))


        edges_possible = choose_random_edges(c1=c1, c2=c2,
                                             m=self.m_high,
                                             rng=self.bm.rng)
        self.edges = edges_possible
    def manages(self, a, b):
        """Return true if two nodes link is managed by this object"""
        if (   (a in self.c1 and b in self.c2)
            or (b in self.c1 and a in self.c2)):
            return True
        return False

    def x_at_t(self, t):
        tau = self.tau
        # Sawtooth profile
        x =  t/float(self.tau) + self.phasefactor
        x = ((x + .5)  % 1) - .5   # x \in (-.5, .5]
        x = abs(x)                 # x \in [0,   .5]
        x *= 2                     # x \in [0,   1]
        ## Cosine profile
        #x = t / float(self.tau)
        #omega = 2*pi*x
        #x = -.5*(math.cos(omega)-1)
        return x
    def m_at_t(self, t):
        x = self.x_at_t(t)
        m = self.m_low + x*(self.m_high-self.m_low)
        debug('Merging: x, m: %s %s %s %s', x, m, self.m_low, self.m_high)
        return m
    def p_at_t(self, t):
        x = self.x_at_t(t)
        p = self.p_low + x*(self.p_high-self.p_low)
        return p

    def t(self, g, t):
        m = self.m_at_t(t)
        edges = self.edges[:int(round(m))]
        #g.add_edges_from(edges)
        for a,b in edges:
            add_edge_nonexists(g, a, b)
    def comms(self, t):
        if self.bm.opts.get('no_det_limit', False):
            x = self.x_at_t(t)
            if x < 1:
                return {self.c_id_1: self.c1,
                        self.c_id_2: self.c2}
        else:
            p = self.p_at_t(t)
            if p < self.p_limit:
                return {self.c_id_1: self.c1,
                        self.c_id_2: self.c2}
        # What is the proper form of this?  We shouldn't reuse c_id_1
        # since it is now a different community.
        return {self.c_id_merged: set.union(self.c1, self.c2)}



class ExpandContract(object):
    def __init__(self, bm, c1, c2, p_in, p_out, tau, fraction=.5,
                 phasefactor=0., c_id_1=0, c_id_2=1):
        self.p_in = p_in
        self.p_out = p_out
        self.c1 = c1
        self.c2 = c2
        self.bm = bm
        self.fraction = fraction
        assert len(c1 & c2) == 0, "Communities must not overlap"
        self.tau = tau
        self.phasefactor = phasefactor
        self.c_id_1 = c_id_1
        self.c_id_2 = c_id_2

        self.order1 = order1 = sorted(shuffled(self.bm.rng, c1))
        self.order2 = order2 = sorted(shuffled(self.bm.rng, c2))

        self.int_1_edges = int_1_edges = { }
        self.ext_2_edges = ext_2_edges = { }
        self.int_2_edges = int_2_edges = { }
        self.ext_1_edges = ext_1_edges = { }
        self.order = order = order1 + list(order2)
        N = len(order)

        c1_minsize = int(round(len(c1)*(1-self.fraction)))
        c1_maxsize = int(round(len(c1)+len(c2)*(self.fraction)))

        with override_numpy_seed(self.bm.rng):
          for i, node in enumerate(order):
            # Internal edges from node to c1
            assert node in order
            if i > 0:
                n_edges = scipy.stats.binom(i, p_in).rvs()
                es = self.bm.rng.sample(order[:i], n_edges)
                int_1_edges[node] = es
                if i >= c1_minsize and n_edges <= 0:
                    raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n1=%s)"%(
                                            self.p_in, i))
                assert node not in es
                debug('ExpandContract: Int. c1=%s, %s %s', self.c_id_1, node, sorted(es))
            else:
                int_1_edges[node] = []
                debug('ExpandContract: Int. c1=%s, %s %s', self.c_id_1, node, [])

            # External edges from node to c1
            if i > 0:
                n_edges = scipy.stats.binom(i, p_out).rvs()
                es = self.bm.rng.sample(order[:i], n_edges)
                ext_1_edges[node] = es
                assert node not in es
                debug('ExpandContract: Ext. c2=%s, %s %s', self.c_id_2, node, sorted(es))

            # Internal edges from node to c2
            if i < N-1:
                n_edges = scipy.stats.binom(N-1-i, p_in).rvs()
                es = self.bm.rng.sample(order[i+1:], n_edges)
                int_2_edges[node] = es
                if i <= c1_maxsize and n_edges <= 0:
                    raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n2=%s)"%(
                                            self.p_in, N-i))
                assert node not in es
                debug('ExpandContract: Int. c2=%s, %s %s', self.c_id_2, node, sorted(es))
            else:
                int_2_edges[node] = []
                debug('ExpandContract: Int. c2=%s, %s %s', self.c_id_2, node, [])


            # External edges from node to c1
            if i < N-1:
                n_edges = scipy.stats.binom(N-1-i, p_out).rvs()
                es = self.bm.rng.sample(order[i+1:], n_edges)
                ext_2_edges[node] = es
                assert node not in es
                debug('ExpandContract: Ext. c1=%s, %s %s', self.c_id_1, node, sorted(es))

        # Check for connectedness of the minimal size subgraphs
        g = nx.Graph()
        for n in order[:c1_minsize]:
            g.add_node(n)
            g.add_edges_from((n, n2) for n2 in int_1_edges[n])
        ccs = nx.connected_components(g)
        print ccs
        if len(ccs) != 1:
            raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n1=%s)"%(
                                    self.p_in, len(g)))
        # c2
        g = nx.Graph()
        for n in order[c1_maxsize:]:
            g.add_node(n)
            g.add_edges_from((n, n2) for n2 in int_2_edges[n])
        ccs = nx.connected_components(g)
        print ccs
        if len(ccs) != 1:
            raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n2=%s)"%(
                                    self.p_in, len(g)))



    def manages(self, a, b):
        """Return true if two nodes link is managed by this object"""
        nodes = set.union(self.c1, self.c2)
        if len(set((a,b)) & nodes) == 2:
            return True
        return False

    def x_at_t(self, t):
        # mod1(x) = x - floor(x)   # gnuplot
        def mod1(x): return x % 1.0

        x = t / float(self.tau) + self.phasefactor
        x = 2 *abs(mod1(x+.5 -.75)-.5)  # .5-1-0-.5 with period 1
        return x
    def c1_size_at_t(self, t):
        x = self.x_at_t(t)
        #bound = int(round(len(self.order)*x))
        low = len(self.c1)*(1-self.fraction)
        high = len(self.c1) + self.fraction*len(self.c2)
        y = x*low + (1-x) * high
        c1 = int(round(y))
        return c1
    def comms(self, t):
        c1 = self.c1_size_at_t(t)
        return {self.c_id_1: self.order[:c1],
                self.c_id_2: self.order[c1:]}

    def t(self, g, t):
        c1size = self.c1_size_at_t(t)
        #print 'merging c1:', x, y, c1
        debug('ExpandContract: t=%s, c1=%s', t, c1size)
        c1 = set(self.order[:c1size])
        c2 = set(self.order[c1size:])

        for i in range(0, c1size):
            n1 = self.order[i]
            #print n1, len(self.int_1_edges[n1]), len(self.ext_2_edges[n1])
            #debug("ExpandContract: adding, c1, Int. %s, %s", self.c_id_1, n1, sorted(self.int_1_edges[n1]))
            for n2 in self.int_1_edges[n1]:
                #print 'a 1 i', n1, n2
                add_edge_nonexists(g, n1, n2)
                #print 'a 2 e', n1, n2
            #debug("ExpandContract: adding, c2, Ext. %s, %s", self.c_id_1, sorted(self.ext_2_edges[n1]))
            for n2 in self.ext_2_edges[n1]:
                if n2 in c2:
                #if n2 > n1:
                    add_edge_nonexists(g, n1, n2)
        for i in range(c1size, len(self.order)):
            n1 = self.order[i]
            #print n1, len(self.ext_1_edges[n1]), len(self.int_2_edges[n1])
            #debug("ExpandContract: adding, c1, Ext. %s, %s", n1, sorted(self.ext_1_edges[n1]))
            #for n2 in self.ext_1_edges[n1]:
            #    if n1 > n2:
            #        add_edge_nonexists(g, n1, n2)
            #debug("ExpandContract: adding, c1=%s, Ext. %s, %s", self.c_id_1, n1, sorted(self.int_2_edges[n1]))
            for n2 in self.int_2_edges[n1]:
                add_edge_nonexists(g, n1, n2)
        # check connectedness (we should never get to this point,
        # should be checked above.  Left for sanity checking, comment
        # this out later).
        if len(nx.connected_components(g.subgraph(self.order[0:c1size]))) != 1:
            ccs = nx.connected_components(g.subgraph(self.order[0:c1size]))
            raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n1=%s, c1=%s, nodes=%s, order=%s, ccs=%s)"%(
                self.p_in, c1size, self.c_id_1, self.order[0:c1size], self.order, ccs))
        if len(nx.connected_components(g.subgraph(self.order[c1size:len(self.order)]))) != 1:
            ccs = nx.connected_components(g.subgraph(self.order[c1size:len(self.order)]))
            raise DisconnectedError("Subgraph is disconnected (ExpandContract, p_in=%f, n2=%s, c2=%s, nodes=%s, order=%s, ccs=%s)"%(
                self.p_in, len(self.order)-c1size, self.c_id_2, self.order[c1size:len(self.order)], self.order, ccs))



re_k = re.compile('k= *([0-9.]+) *')
class _StdBase(Benchmark):
    _default_p_in = 'k=16'
    _default_p_out = 'k=0'
    def __init__(self, p_in=1., p_out=0., n=32, q=4, tau=100,
                 opts={}):
        self.rng = random.Random(opts.get('seed', None))
        self.opts = opts

        if isinstance(p_in, str) and re_k.match(p_in):
            p_in  = float(re_k.match(p_in).group(1)) / (n-1)
        if isinstance(p_out, str) and re_k.match(p_out):
            p_out = float(re_k.match(p_out).group(1)) / n
        self.p_in = p_in
        self.p_out = p_out



class StdMerge(_StdBase):
    def __init__(self, p_in=_StdBase._default_p_in, p_out=_StdBase._default_p_out,
                 n=32, q=4, tau=100, opts={}):
        super(StdMerge, self).__init__(p_in=p_in, p_out=p_out, n=n, q=q,
                                       tau=tau, opts=opts)
        p_in = self.p_in
        p_out = self.p_out

        if q%2 != 0:
            raise ValueError("q must be a multiple of two (given: q=%s)"%q)

        cs = [set(range(n*i, n*(i+1))) for i in range(q)]

        managers = [ ]
        for i in range(q//2):
            c0 = 2*i
            c1 = 2*i+1
            managers.append(
                Merging(self, cs[c0], cs[c1],
                        p_high=p_in, p_low=p_out, tau=tau,
                        phasefactor=i/float(q//2),
                        c_id_1=c0, c_id_2=c1, c_id_merged=q+i))
            managers.append(Static(self, cs[c0], p=p_in))
            managers.append(Static(self, cs[c1], p=p_in))
            for j in range(i+1, q//2):
                d0 = 2*j
                d1 = 2*j+1
                managers.append(
                    Static(self, c1=cs[c0]|cs[c1],
                                 c2=cs[d0]|cs[d1], p=p_out))
        self.managers = managers
        self.g = g = nx.Graph()

        for c in cs:
            for n in c:
                g.add_node(n)

class StdGrow(_StdBase):
    def __init__(self, p_in=_StdBase._default_p_in, p_out=_StdBase._default_p_out,
                 n=32, q=4, tau=100, opts={}):
        super(StdGrow, self).__init__(p_in=p_in, p_out=p_out, n=n, q=q,
                                       tau=tau, opts=opts)
        p_in = self.p_in
        p_out = self.p_out

        if q%2 != 0:
            raise ValueError("q must be a multiple of two (given: q=%s)"%q)

        cs = [set(range(n*i, n*(i+1))) for i in range(q)]

        managers = [ ]
        for i in range(q//2):
            c0 = 2*i
            c1 = 2*i+1
            managers.append(
                ExpandContract(self, cs[c0], cs[c1],
                               p_in=p_in, p_out=p_out, tau=tau,
                               phasefactor=i/float(q//2),
                               c_id_1=c0, c_id_2=c1))
            for j in range(i+1, q//2):
                d0 = 2*j
                d1 = 2*j+1
                managers.append(
                    Static(self, cs[c0]|cs[c1],
                                 cs[d0]|cs[d1], p=p_out))
        self.managers = managers

        self.g = g = nx.Graph()
        # Add all initial nodes to the graph
        for c in cs:
            for n in c:
                g.add_node(n)

class StdMixed(_StdBase):
    def __init__(self, p_in=_StdBase._default_p_in, p_out=_StdBase._default_p_out,
                 n=32, q=4, tau=100, opts={}):
        super(StdMixed, self).__init__(p_in=p_in, p_out=p_out, n=n, q=q,
                                       tau=tau, opts=opts)
        p_in = self.p_in
        p_out = self.p_out

        if q%4 != 0:
            raise ValueError("q must be a multiple of four (given: q=%s)"%q)

        cs = [set(range(n*i, n*(i+1))) for i in range(q)]

        managers = [ ]
        for i in range(q//4):
            c0, c1, c2, c3 = 4*i, 4*i+1, 4*i+2, 4*i+3
            managers.append(
                Merging(self, cs[c0], cs[c1],
                        p_high=p_in, p_low=p_out, tau=tau,
                        phasefactor=i/float(q//4),
                        c_id_1=c0, c_id_2=c1, c_id_merged=q+i))
            managers.append(Static(self, cs[c0], p=p_in))
            managers.append(Static(self, cs[c1], p=p_in))
            managers.append(
                ExpandContract(self, cs[c2], cs[c3],
                               p_in=p_in, p_out=p_out, tau=tau,
                               phasefactor=i/float(q//4),
                               c_id_1=c2, c_id_2=c3))
            managers.append(Static(self, cs[c0]|cs[c1],
                                         cs[c2]|cs[c3], p=p_out))
            for j in range(i+1, q//4):
                d0, d1, d2, d3 = 4*j, 4*j+1, 4*j+2, 4*j+3
                managers.append(
                    Static(self, cs[c0]|cs[c1]|cs[c2]|cs[c3],
                                 cs[d0]|cs[d1]|cs[d2]|cs[d3], p=p_out),
                    )
        self.managers = managers

        self.g = g = nx.Graph()
        # Add all initial nodes to the graph
        for c in cs:
            for n in c:
                g.add_node(n)


def main_argv(argv=sys.argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("bm_model", help="benchmark model to simulate")

    parser.add_argument("output", help="Output prefix", nargs='?')
    parser.add_argument("--t",  type=int, help="Maximum time to simulate",
                        default=100)

    group =parser.add_argument_group(title="Model parameters",description=None)
    parser.add_argument("--q", help="Number of communities", type=int)
    parser.add_argument("--n", help="", type=int)
    parser.add_argument("--p_in", help="Internal edge density", type=float)
    parser.add_argument("--p_out", help="External edge density",type=float)
    parser.add_argument("--k_in",  type=float)
    parser.add_argument("--k_out", type=float)
    parser.add_argument("--tau",   type=int)
    parser.add_argument("--graph-format", default='edgelist', help="How to write graph, choices='edgelist', 'pajek'.")
    parser.add_argument("--comm-format", default='bynode', help="How to write communities, choices='oneline', 'bynode', 'pajek'.")
    #parser.add_argument("--", help="", type=int, default=)
    parser.add_argument("--seed",  default=None, help="Random seed")
    parser.add_argument("--no-det-limit",  action='store_true',
                        help="No detectability limit")
    parser.add_argument("--Gnm",  action='store_true',
                        help="Use Gnm random graph ensemble instead of Gnp.  Only works for merging, not grow/shrink.")
    model_params_names = ['q', 'n', 'p_in', 'p_out', 'tau', ]

    print argv
    args = parser.parse_args(args=argv[1:])

    model_params = dict((name, getattr(args, name))
                        for name in model_params_names
                        if getattr(args, name) is not None)
    print model_params
    if args.k_in is not None:
        model_params['p_in']  = 'k=%f'%args.k_in
    if args.k_out is not None:
        model_params['p_out'] = 'k=%f'%args.k_out

    return (get_model(args.bm_model, opts=args.__dict__, **model_params),
            args)

def get_model(name=None, **kwargs):
    """Return a given model names, instantiated with **kwargs"""
    bm = globals()[name]
    bm = bm(**kwargs)
    return bm

def main(argv=sys.argv):
    """Main entry point from command line."""
    bm, args = main_argv(argv)
    run(bm, maxt=args.t, output=args.output,
        graph_format=args.graph_format,
        comm_format=args.comm_format,
        opts=args.__dict__)
    return bm

def run(bm, maxt=100, output=None, graph_format='edgelist',
        comm_format='bynode',
        opts={}):
    """Main loop to do a running."""
    for t in range(maxt+1):
        g = bm.t(t)
        comms = bm.comms(t)
        if output:
            prefix = output + '.t%05d'%t
            # write graph
            if graph_format == 'edgelist':
                nx.write_edgelist(g, prefix+'.graph', data=False)
            elif graph_format == 'pajek':
                for n in g:
                    g.node[n]['id'] = n+1
                nx.write_pajek(g, prefix+'.graph')
            elif graph_format == 'null':
                pass
            else:
                try:
                    graphwriter = getattr(nx, 'write_'+graph_format)
                except AttributeError:
                    raise ValueError("Unknown graph format: %s"%graph_format)
                graphwriter(g, prefix+'.graph')
            # write communities, in desired format.
            if comm_format == 'oneline':  comm_writer = write_comms_oneline
            elif comm_format == 'bynode': comm_writer = write_comms_bynode
            elif comm_format == 'pajek':  comm_writer = write_comms_pajek
            elif comm_format == 'null':   comm_writer = None
            else:
                raise ValueError("Unknown comm format: %s"%comm_format)
            # Delay opening the file until here, so we can skip it
            # using the 'null' option.
            if comm_writer:
                f = open(prefix+'.comms', 'w')
                label = 't=%s, command line: %s'%(t, ' '.join(sys.argv))
                comm_writer(f, comms, label)
        print t, len(g), g.number_of_edges(), len(comms), sorted(comms.keys())


def write_comms_oneline(f, comms, label=None):
    """Write communities, one line per community."""
    if label:
        print >> f, '#', label.replace('\n', ' ')
    print >> f, '#', time.ctime()
    print >> f, '# Format: "node_id node_id node_id ...", one line per community.'
    for cname, cnodes in comms.iteritems():
        print >> f, "# label: %s"%cname
        print >> f, " ".join(str(x) for x in cnodes)
def write_comms_bynode(f, comms, label=None):
    """Write communities, lines with 'node comm' pairs"""
    if label:
        print >> f, '#', label.replace('\n', ' ')
    print >> f, '#', time.ctime()
    print >> f, "# Format: node_id cmty_id"
    for cname, cnodes in comms.iteritems():
        for node in cnodes:
            print >> f, node, cname

def write_comms_pajek(f, comms, label=None):
    """Write communities, lines with 'node comm' pairs"""
    if label:
        print >> f, '#', label.replace('\n', ' ')
    print >> f, '#', time.ctime()
    print >> f, "# Format: cmty_id in node order"
    print >> f, "*vertices"
    nodecmtys = { }
    for cname, cnodes in comms.iteritems():
        for node in cnodes:
            nodecmtys[node] = cname
    for node in sorted(nodecmtys):
        print >> f, nodecmtys[node]


if __name__ == "__main__":
    main()
