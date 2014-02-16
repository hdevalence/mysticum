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

if __name__ == '__main__':
    unittest.main()
