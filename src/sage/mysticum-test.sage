#!/usr/bin/env sage
"""
Test suite for mysticum.sage

(C) 2014 David Wehlau / Royal Military College of Canada

Author: Henry de Valence
"""

import mysticum as myst
import unittest

class TestMysticum(unittest.TestCase):
    def setUp(self):
        # Points on the unit circle
        a = vector([5/13, 12/13, 1])
        b = vector([1, 0, 1])
        c = vector([8/17, -15/17, 1])
        d = vector([-4/5, -3/5, 1])
        e = vector([-40/41, -9/41, 1])
        f = vector([-3/5, 4/5, 1])
        self.hexagon = (a,b,c,d,e,f)

    def test_line_through_two_points(self):
        p1 = vector(QQ,[1,0,1])
        p2 = vector(QQ,[0,1,1])
        l = myst.line_through_points(p1,p2)
        self.assertEqual(l, vector(QQ,[1,1,-1]))

    def test_line_through_one_point(self):
        p1 = [1,0,0]
        with self.assertRaises(myst.UnderDeterminedError):
            myst.line_through_points(p1)

    def test_line_through_three_points(self):
        p1 = [1,0,1]
        p2 = [0,1,1]
        p3 = [1,1,1]
        with self.assertRaises(myst.OverDeterminedError):
            myst.line_through_points(p1,p2,p3)

    def test_line_intersection(self):
        l1 = [1,1,-1] # x + y = 1
        l2 = [1,1/2,-1] # x + 1/2y = 1
        p = myst.line_intersection(l1,l2)
        self.assertEqual(p, vector(QQ,[1,0,1]))

    def test_pascal_line(self):
        l = myst.pascal_line(*self.hexagon)
        self.assertEqual(l, vector(QQ, [1, 18/179, 847/179]))

    def test_pascal_lines(self):
        G = SymmetricGroup(6)
        g = G((1,3,5,2,6,4))
        lines = myst.pascal_lines(*self.hexagon)
        self.assertEqual(lines[g], vector(QQ, (1, -221/115, 17/115)))

    def test_60_hexagons(self):
        G = SymmetricGroup(6)
        hexes = myst.conj_class_up_to_inverses(G((1,2,3,4,5,6)))
        self.assertEqual(len(hexes), 60)

    def test_cycle_edges(self):
        G = SymmetricGroup(6)
        e = myst.cycle_edges(G((1,2,3,4,5,6)))
        self.assertEqual(e, set([(1,2),(5,6),(4,5),(2,3),(1,6),(3,4)]))

    def test_cycles_disjoint(self):
        # N.B. There are some issues to be resolved about
        # this test, since we only care about things up to inverses.
        G = SymmetricGroup(6)
        disjoints = myst.cycles_disjoint(G((1,2,3,4,5,6)))
        self.assertEqual(disjoints, set([G((1,5,3,6,2,4)),
                                         G((1,3,5,2,6,4)),
                                         G((1,3,6,4,2,5))]))
    def test_3_disjoint(self):
        G = SymmetricGroup(6)
        sixcycles = G.conjugacy_class(G((1,2,3,4,5,6)))
        for g in sixcycles:
            self.assertEqual(len(myst.cycles_disjoint(g)),3)

    def test_apply_cycle(self):
        G = SymmetricGroup(6)
        g = G((1,3,5,2,6,4))
        a,b,c,d,e,f = self.hexagon
        self.assertEqual(myst.apply_cycle(g, self.hexagon), (a,c,e,b,f,d))

    def test_kirkman_node(self):
        G = SymmetricGroup(6)
        g = G((1,2,3,4,5,6))
        node = myst.kirkman_node(g,self.hexagon)
        self.assertEqual(node, vector(QQ,(1, 38/153, -541/153)))

    def test_60_kirkman_nodes(self):
        nodes = myst.kirkman_nodes(self.hexagon)
        self.assertEqual(len(nodes), 60)

if __name__ == '__main__':
    unittest.main()
