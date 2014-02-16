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
    hexagons = set()
    for g in G.conjugacy_class(G((1,2,3,4,5,6))):
        if g^(-1) not in hexagons:
            hexagons.add(g)
    lines = dict()
    for g in hexagons:
        lines[g] = pascal_line(*g(points))
    return lines

