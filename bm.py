import math
import random
import sys

import networkx as nx
import scipy.stats

import logging
logger = logging.getLogger(__name__)
debug = logger.debug
info = logger.info

class DisconnectedError(Exception):
    pass


class Benchmark(object):
    def __init__(self, p_in=1, p_out=0, tau=100):
        self.rng = random.Random()

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

    def t(self, t):
        g = self.g.copy()
        for mgr in self.managers:
            mgr.t(g, t)
        return g
    def comms(self, t):
        comms = [ ]
        for mrg in self.managers:
            comms.extend(mrg.comms(t))
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
                n1 = random.sample(c1, 1)[0]#sets support sample but not choice
                n2 = random.sample(c2, 1)[0]
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

        # Number of total possible edges
        if c2 is None:
            n_links = len(c1) * (len(c1)-1) / 2
        else:
            n_links = len(c1) * len(c2)
        # Number of actual edges
        n_edges = scipy.stats.binom(n_links, p).rvs()

        self.edges_active = choose_random_edges(c1=c1, c2=c2, m=n_edges,
                                           rng=self.bm.rng)

    def t(self, g, t):
        for a,b in self.edges_active:
            add_edge_nonexists(g, a, b)
        #g.add_edges_from(self.edges_active)

        if self.c2 is None:
            if len(nx.connected_components(g.subgraph(self.c1))) != 1:
                raise DisconnectedError("Subgraph is disconnected (Static object)")
    def comms(self, t):
        return [ ]
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
    def __init__(self, bm, c1, c2, p_low, p_high, tau, phasefactor=0.):
        self.bm = bm
        self.n_links = n_links = len(c1) * len(c2)
        self.c1 = c1
        self.c2 = c2

        self.p_low = p_low
        self.p_high = p_high
        self.m_low  = scipy.stats.binom(n_links, p_low ).rvs()
        self.m_high = scipy.stats.binom(n_links, p_high).rvs()
        self.tau = tau
        self.phasefactor = phasefactor

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
        return m


    def t(self, g, t):
        m = self.m_at_t(t)
        edges = self.edges[:int(round(m))]
        #g.add_edges_from(edges)
        for a,b in edges:
            add_edge_nonexists(g, a, b)
    def comms(self, t):
        x = self.x_at_t(t)
        if x < 1:
            return [self.c1, self.c2]
        return [set.union(self.c1, self.c2)]



class ExpandContract(object):
    def __init__(self, bm, c1, c2, p_in, p_out, tau, fraction=.5,
                 phasefactor=0.):
        self.c1 = c1
        self.c2 = c2
        self.bm = bm
        self.fraction = fraction
        assert len(c1 & c2) == 0, "Communities must not overlap"
        self.tau = tau
        self.phasefactor = phasefactor

        self.order1 = order1 = sorted(shuffled(self.bm.rng, c1))
        self.order2 = order2 = sorted(shuffled(self.bm.rng, c2))

        self.int_1_edges = int_1_edges = { }
        self.ext_2_edges = ext_2_edges = { }
        self.int_2_edges = int_2_edges = { }
        self.ext_1_edges = ext_1_edges = { }
        self.order = order = order1 + list(reversed(order2))
        N = len(order)

        for i, node in enumerate(order):
            # Internal edges from node to c1
            assert node in order
            if i > 0:
                n_edges = scipy.stats.binom(i, p_in).rvs()
                es = self.bm.rng.sample(order[:i], n_edges)
                int_1_edges[node] = es
                #if i >= self.fraction*len(c1): assert n_edges > 0
                assert node not in es
            else:
                int_1_edges[node] = []

            # External edges from node to c1
            if i > 0:
                n_edges = scipy.stats.binom(i, p_out).rvs()
                es = self.bm.rng.sample(order[:i], n_edges)
                ext_1_edges[node] = es
                assert node not in es

            # Internal edges from node to c2
            if i < N-1:
                n_edges = scipy.stats.binom(N-1-i, p_in).rvs()
                es = self.bm.rng.sample(order[i+1:], n_edges)
                int_2_edges[node] = es
                #if i <= N-self.fraction*len(c2): assert n_edges > 0
                assert node not in es
            else:
                int_2_edges[node] = []


            # External edges from node to c2
            if i < N-1:
                n_edges = scipy.stats.binom(N-1-i, p_out).rvs()
                es = self.bm.rng.sample(order[i+1:], n_edges)
                ext_2_edges[node] = es
                assert node not in es
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
        return [self.order[:c1], self.order[c1:]]

    def t(self, g, t):
        c1size = self.c1_size_at_t(t)
        #print 'merging c1:', x, y, c1
        c1 = set(self.order[:c1size])
        c2 = set(self.order[c1size:])

        for i in range(0, c1size):
            n1 = self.order[i]
            #print n1, len(self.int_1_edges[n1]), len(self.ext_2_edges[n1])
            for n2 in self.int_1_edges[n1]:
                #print 'a 1 i', n1, n2
                add_edge_nonexists(g, n1, n2)
                #print 'a 2 e', n1, n2
            for n2 in self.ext_2_edges[n1]:
                if n2 in c2:
                #if n2 > n1:
                    add_edge_nonexists(g, n1, n2)
        for i in range(c1size, len(self.order)):
            n1 = self.order[i]
            #print n1, len(self.ext_1_edges[n1]), len(self.int_2_edges[n1])
            #for n2 in self.ext_1_edges[n1]:
            #    if n1 > n2:
            #        add_edge_nonexists(g, n1, n2)
            for n2 in self.int_2_edges[n1]:
                add_edge_nonexists(g, n1, n2)


class StdMerge(Benchmark):
    def __init__(self, p_in=1., p_out=0., n=32, q=4, tau=100):
        self.rng = random.Random()

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
                        phasefactor=i/float(q//2)))
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

class StdGrow(Benchmark):
    def __init__(self, p_in=1, p_out=0, n=32, q=4, tau=100):
        self.rng = random.Random()

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
                               phasefactor=i/float(q//2)))
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

class StdMixed(Benchmark):
    def __init__(self, p_in=1, p_out=0, n=32, q=4, tau=100):
        self.rng = random.Random()

        if q%4 != 0:
            raise ValueError("q must be a multiple of four (given: q=%s)"%q)

        cs = [set(range(n*i, n*(i+1))) for i in range(q)]

        managers = [ ]
        for i in range(q//4):
            c0, c1, c2, c3 = 4*i, 4*i+1, 4*i+2, 4*i+3
            managers.append(
                Merging(self, cs[c0], cs[c1],
                        p_high=p_in, p_low=p_out, tau=tau,
                        phasefactor=i/float(q//4)))
            managers.append(Static(self, cs[c0], p=p_in))
            managers.append(Static(self, cs[c1], p=p_in))
            managers.append(
                ExpandContract(self, cs[c2], cs[c3],
                               p_in=p_in, p_out=p_out, tau=tau,
                               phasefactor=i/float(q//4)))
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
    parser.add_argument("--k_in",  type=int)
    parser.add_argument("--k_out", type=int)
    parser.add_argument("--tau",   type=int)
    #parser.add_argument("--", help="", type=int, default=)
    model_params_names = ['q', 'n', 'p_in', 'p_out', 'tau']

    print argv
    args = parser.parse_args(args=argv[1:])

    model_params = dict((name, getattr(args, name))
                        for name in model_params_names
                        if getattr(args, name) is not None)
    print model_params
    if args.k_in is not None:
        model_params['p_in']  = args.k_in  / float(args.n - 1)
    if args.k_out is not None:
        model_params['p_out'] = args.k_out / float(args.n)

    return get_model(args.bm_model, **model_params)

def get_model(name=None, **kwargs):
    """Return a given model names, instantiated with **kwargs"""
    bm = globals()[name]
    bm = bm(**kwargs)
    return bm

def main(argv=sys.argv):
    bm = main_argv(argv)
    for t in range(args.t + 1):
        g = bm.t(t)
        if args.output:
            prefix = args.output + '.t%05d'%t
            nx.write_edgelist(g, prefix+'.edges', data=False)
            comms = bm.comms(t)
            f = open(prefix+'.comms', 'w')
            write_comms(f, comms)

        print t, len(g), g.number_of_edges()


def write_comms(f, comms):
    for nodes in comms:
        print >> f, " ".join(str(x) for x in nodes)


if __name__ == "__main__":
    main()
