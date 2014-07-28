import math
import random

import networkx as nx
import scipy.stats

class Benchmark(object):
    def __init__(self, p_low=0.1, p_high=.5, tau=50):
        self.rng = random.Random()

        n = 32

        self.c1 = c1 = set(range(0, n  ))
        self.c1 = c2 = set(range(n, n*2))

        nodes = c1 | c2

        managers = [Static(self, c1, p=p_high),
                    Static(self, c2, p=p_high),
                    Merging(self, c1, c2, p_low=p_low, p_high=p_high, tau=tau),
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
        g.add_edges_from(self.edges_active)

        if self.c2 is None:
            assert len(nx.connected_components(g.subgraph(self.c1))) == 1


class Merging(object):
    def __init__(self, bm, c1, c2, p_low, p_high, tau):
        self.bm = bm
        self.n_links = n_links = len(c1) * len(c2)

        self.p_low = p_low
        self.p_high = p_high
        self.m_low  = scipy.stats.binom(n_links, p_low ).rvs()
        self.m_high = scipy.stats.binom(n_links, p_high).rvs()
        self.tau = tau

        edges_possible = choose_random_edges(c1=c1, c2=c2,
                                             m=self.m_high,
                                             rng=self.bm.rng)
        self.edges = edges_possible

    def m_at_t(self, t):
        tau = self.tau
        # Sawtooth profile
        x =  t/float(self.tau)
        x = ((x + .5)  % 1) - .5
        x = abs(x)
        x *= 2  # scale to max of 1 / min of 2
        m = self.m_low + x*(self.m_high-self.m_low)
        return m
        # Cosine profile
        x = t / float(self.tau)
        omega = 2*pi*x
        x = -.5*(math.cos(omega)-1)
        m = self.m_low + x*(self.m_high-self.m_low)
        return m


    def t(self, g, t):
        m = self.m_at_t(t)
        edges = self.edges[:int(round(m))]
        g.add_edges_from(edges)



if __name__ == "__main__":
    bm = Benchmark()
    for t in range(50):
        g = bm.t(t)
        print t, len(g), g.number_of_edges()

