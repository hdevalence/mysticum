#!/usr/bin/env sage
# coding: utf-8
"""
Code to compute the complete Mysticum Hexagrammaticum.

(C) 2014 David Wehlau / Royal Military College of Canada

Author: Henry de Valence
"""

from instance_memoize import memoize

class OverDeterminedError(Exception):
    """Error raised by over-determined constraints."""
    pass

class UnderDeterminedError(Exception):
    """Error raised by under-determined constraints."""
    pass

def line_through_points(*points):
    """
    Given points in projective space, compute the line through them.
    """
    B = column_matrix(points).kernel().basis()
    if len(B) == 1:
        return B[0]
    elif len(B) == 0:
        raise OverDeterminedError
    else:
        raise UnderDeterminedError

def line_intersection(*lines):
    """
    Given lines in projective space, compute their point of intersection.
    """
    return line_through_points(*lines)

def conj_class_up_to_inverses(g):
    """
    Return the conjugacy class of g, with either h or h^-1
    discarded for every h in the conjugacy class.
    """
    conj_class = set()
    for h in g.parent().conjugacy_class(g):
        if h^(-1) not in conj_class:
            conj_class.add(h)
    return conj_class

def cycle_edges(g):
    """
    Given g, an n-cycle, compute the set of edges of g.
    If g = (abcdef), the edges are ab. bc, ..., fa.
    """
    e = set()
    x = 1
    for _ in g.orbit(x):
        # Edges are unordered, so impose an ordering to
        # ensure a unique representation.
        if x < g(x):
            e.add((x,g(x)))
        else:
            e.add((g(x),x))
        x = g(x)
    return e

def cycles_disjoint(g):
    """
    Given g, an n-cycle, return the cycles disjoint from g,
    i.e., those with no common edges.
    """
    e = cycle_edges(g)
    return set((h for h in conj_class_up_to_inverses(g)
                if cycle_edges(h).isdisjoint(e)))

def cycles_squaring_to(s):
    """
    Given s, a cycle of shape (3,3), return the 6-cycles
    (up to inverses) whose squares are s.
    """
    G = SymmetricGroup(6)
    g = G((1,2,3,4,5,6))
    sixcycles = G.conjugacy_class(g)
    return set((h for h in sixcycles if h^2 == s))

def apply_cycle(g,hex):
    """
    Reorders the points in the hexagon hex according to the cycle g,
    and returns the new hexagon.
    """
    new_hex = [None]*6
    x = 1
    for i in range(6):
        new_hex[i] = hex[x-1]
        x = g(x)
    return tuple(new_hex)

def cycles_commuting_with(C, g):
    """
    Return cycles (up to inverses) in C which commute with g.
    """
    return (h for h in C if (g*h == h*g) and h < h.inverse())

class Mysticum:
    """
    Object representing the 95 lines and 95 points associated to
    a given hexagon.
    """
    G = SymmetricGroup(6)

    def __init__(self, a, b, c, d, e, f):
        """
        Constructs the Mysticum through the given points,
        which are assumed to lie on a conic in PP^2. The points
        are given as 3-vectors in homogeneous coordinates.
        """
        self.hexagon = (a,b,c,d,e,f)

    @memoize
    def pascal_line(self, g):
        """
        Get the pascal line L(g), for g a 6-cycle.
        """
        hex = apply_cycle(g, self.hexagon)
        p0 = line_intersection(line_through_points(hex[0],hex[1]),
                               line_through_points(hex[3],hex[4]))
        p1 = line_intersection(line_through_points(hex[1],hex[2]),
                               line_through_points(hex[4],hex[5]))
        p2 = line_intersection(line_through_points(hex[2],hex[3]),
                               line_through_points(hex[5],hex[0]))
        return line_through_points(p0,p1,p2)

    @memoize
    def kirkman_node(self, g):
        """
        Compute the Kirkman node N(g), for g a 6-cycle.
        """
        lines = [self.pascal_line(h) for h in cycles_disjoint(g)]
        return line_intersection(*lines)

    @memoize
    def steiner_node(self, g):
        """
        Compute the Steiner node N(g) for g of shape (3,3).
        """
        lines = [self.pascal_line(h) for h in cycles_squaring_to(g)]
        return line_intersection(*lines)

    @memoize
    def cayley_line(self, g):
        """
        Compute the Cayley line L(g) for g of shape (3,3).
        """
        nodes = [self.kirkman_node(h) for h in cycles_squaring_to(g)]
        return line_through_points(*nodes)

    @memoize
    def pluecker_line(self, g):
        """
        Compute the Plücker line L(g) for g of shape (2,2,2).
        """
        C = self.G("(1,2,3)(4,5,6)").conjugacy_class()
        nodes = [self.steiner_node(h) for h in cycles_commuting_with(C, g)]
        return line_through_points(*nodes)

