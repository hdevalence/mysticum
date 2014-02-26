#!/usr/bin/env sage
"""
Code to compute the complete Mysticum Hexagrammaticum.

(C) 2014 David Wehlau / Royal Military College of Canada

Author: Henry de Valence
"""

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

def pascal_line(*hexagon):
    """
    Given a 6-tuple of points, compute the Pascal line through these points.
    """
    p0 = line_intersection(line_through_points(hexagon[0],hexagon[1]),
                           line_through_points(hexagon[3],hexagon[4]))
    p1 = line_intersection(line_through_points(hexagon[1],hexagon[2]),
                           line_through_points(hexagon[4],hexagon[5]))
    p2 = line_intersection(line_through_points(hexagon[2],hexagon[3]),
                           line_through_points(hexagon[5],hexagon[0]))
    return line_through_points(p0,p1,p2)

def pascal_lines(*points):
    """
    Given a 6-tuple of points, compute all the Pascal lines through these points.
    """
    G = SymmetricGroup(6)
    # Hexagons <--> 6-cycles, up to inverses
    hexagons = conj_class_up_to_inverses(G((1,2,3,4,5,6)))
    lines = dict()
    for g in hexagons:
        lines[g] = pascal_line(*apply_cycle(g,points))
    return lines

def kirkman_node(g, points):
    """
    Compute the Kirkman node for points, with ordering given by g.

    The Kirkman nodes N(g) are indexed by 6-cycles; they are
    the intersection of the Pascal lines corresponding to the
    6-cycles disjoint from a given 6-cycle.
    """
    lines = [pascal_line(*apply_cycle(h,points)) for h in cycles_disjoint(g)]
    return line_intersection(*lines)

def kirkman_nodes(points):
    """
    Compute all Kirkman nodes for points, indexed by 6-cycles.
    """
    G = SymmetricGroup(6)
    cycles = conj_class_up_to_inverses(G((1,2,3,4,5,6)))
    nodes = dict()
    for g in cycles:
        nodes[g] = kirkman_node(g,points)
    return nodes

def steiner_node(g, points):
    """
    Compute the Steiner node corresponding to g, a permutation
    of shape (3,3).
    """
    lines = [pascal_line(*apply_cycle(h,points))
             for h in cycles_squaring_to(g)]
    return line_intersection(*lines)

def steiner_nodes(points):
    """
    Compute all Steiner nodes for points, indexed by (3,3)-cycles.
    """
    G = SymmetricGroup(6)
    cycles = conj_class_up_to_inverses(G("(1,2,3)(4,5,6)"))
    nodes = dict()
    for g in cycles:
        nodes[g] = steiner_node(g,points)
    return nodes

def cayley_line(g, points):
    """
    Compute the Cayley line corresponding to the (3,3)-cycle g.
    """
    nodes = [kirkman_node(h, points) for h in cycles_squaring_to(g)]
    return line_through_points(*nodes)

def cayley_lines(points):
    """
    Compute the Cayley lines for points, indexed by (3,3)-cycles.
    """
    G = SymmetricGroup(6)
    cycles = conj_class_up_to_inverses(G("(1,2,3)(4,5,6)"))
    lines = dict()
    for g in cycles:
        lines[g] = cayley_line(g,points)
    return lines

