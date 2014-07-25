import math
import random

import networkx as nx
import scipy.stats

class Benchmark(object):
    def __init__(self, p_low=0.1, p_high=.5, tau=50):
        self.rng = random.Random()


        c1 = set(range( 0,32))
        c2 = set(range(32,64))

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
            mgr.t(g)
        return g


class Static(object):
    def __init__(self, bm, c1, c2=None, p=None):
        self.bm = bm
        assert p is not None
        if c2 is not None:
            edges = set(frozenset((a, b))  for a in c1  for b in c2)
            assert len(edges) == len(c1) * len(c2), "overlap between c1 and c2?"
        else:
            lst = list(c1)
            edges = set(frozenset((lst[i], lst[j]))  for i in range(len(lst)) for j in range(i+1, len(lst)) )
            assert len(edges) == len(c1) * (len(c1)-1) / 2

        edges = list(edges)
        self.bm.rng.shuffle(edges)

        n_edges = scipy.stats.binom(len(edges), p).rvs()

        self.edges_active = edges[:n_edges]

    def t(self, g):
        g.add_edges_from(self.edges_active)


class Merging(object):
    def __init__(self, bm, c1, c2, p_low, p_high, tau):
        self.bm = bm
        pairs = len(c1) * len(c2)

        self.p_low = p_low
        self.p_high = p_high
        self.tau = tau


        edges_possible = set(frozenset((a, b))  for a in c1  for b in c2)
        assert len(edges_possible) == len(c1) * len(c2), "overlap between c1 and c2?"
        edges_possible = list(edges_possible)
        self.bm.rng.shuffle(edges_possible)

        self.edges = edges_possible

    def density_at_t(self, t):
        tau = self.tau
        x =  t/float(self.tau)   #
        x = ((x + .5)  % 1) - .5
        return abs(x)
        


        
        return self.p_low + (self.p_high-self.p_low) * -.5*(math.cos(omega)-1)


    def t(self, g):
        rho = self.density_at_t(t)

        n_edges = scipy.stats.binom(len(self.edges), rho).rvs()
        n_edges = int(rho*len(self.edges))
        edges = self.edges[:n_edges]

        g.add_edges_from(edges)



if __name__ == "__main__":
    bm = Benchmark()
    for t in range(50):
        g = bm.t(t)
        print len(g), g.number_of_edges()

