import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestCells(unittest.TestCase):

    def test_parent(self):

        region = openmoc.Union()
        self.assertEqual(region.getParentRegion(), None)
        region.setParentRegion(openmoc.Union())
        self.assertNotEqual(region.getParentRegion(), None)

    def test_union(self):

        union = openmoc.Union()

        # Define surfaces
        p1x = openmoc.XPlane(x = 3)
        p2x = openmoc.XPlane(x = -2)
        p1y = openmoc.YPlane(y = 1)
        p2y = openmoc.YPlane(y = -0.5)
        p1z = openmoc.ZPlane(z = 8)
        p2z = openmoc.ZPlane(z = -4)

        # Define boundary types
        p1x.setBoundaryType(openmoc.VACUUM)
        p2x.setBoundaryType(openmoc.REFLECTIVE)
        p1y.setBoundaryType(openmoc.PERIODIC)
        p2y.setBoundaryType(openmoc.REFLECTIVE)
        p1z.setBoundaryType(openmoc.PERIODIC)
        p2z.setBoundaryType(openmoc.VACUUM)

        # Define halfspaces
        h1x = openmoc.Halfspace(-1, p1x)
        h2x = openmoc.Halfspace(+1, p2x)
        h1y = openmoc.Halfspace(-1, p1y)
        h2y = openmoc.Halfspace(+1, p2y)
        h1z = openmoc.Halfspace(-1, p1z)
        h2z = openmoc.Halfspace(+1, p2z)

        # Test getMaxXYZ + MaxXYZboundary
        union.addNode(h1x, True)
        self.assertEqual(union.getMaxX(), 3)
        self.assertEqual(union.getMaxXBoundaryType(), openmoc.VACUUM)
        #union.removeHalfspace(p1x, -1)  #FIXME seg faults

        union = openmoc.Union()
        union.addNode(h1y)
        self.assertEqual(union.getMaxY(), 1)
        self.assertEqual(union.getMaxYBoundaryType(), openmoc.PERIODIC)
        #union.removeHalfspace(p1y, -1)

        union = openmoc.Union()
        union.addNode(h1z)
        self.assertEqual(union.getMaxZ(), 8)
        self.assertEqual(union.getMaxZBoundaryType(), openmoc.PERIODIC)
        #union.removeHalfspace(p1z, -1)

        # Test getMinXYZ + MinXYZboundary
        union = openmoc.Union()
        union.addNode(h2x)
        self.assertEqual(union.getMinX(), -2)
        self.assertEqual(union.getMinXBoundaryType(), openmoc.REFLECTIVE)
        #union.removeHalfspace(h2x)

        union = openmoc.Union()
        union.addNode(h2y)
        self.assertEqual(union.getMinY(), -0.5)
        self.assertEqual(union.getMinYBoundaryType(), openmoc.REFLECTIVE)
        #union.removeHalfspace(h2y)

        union = openmoc.Union()
        union.addNode(h2z)
        self.assertEqual(union.getMinZ(), -4)
        self.assertEqual(union.getMinZBoundaryType(), openmoc.VACUUM)
        #union.removeHalfspace(h2z)

    def test_intersection(self):

        intersection = openmoc.Intersection()

        # Define surfaces
        p1x = openmoc.XPlane(x = 3)
        p2x = openmoc.XPlane(x = -2)
        p1y = openmoc.YPlane(y = 1)
        p2y = openmoc.YPlane(y = -0.5)
        p1z = openmoc.ZPlane(z = 8)
        p2z = openmoc.ZPlane(z = -4)

        # Define boundary types
        p1x.setBoundaryType(openmoc.VACUUM)
        p2x.setBoundaryType(openmoc.REFLECTIVE)
        p1y.setBoundaryType(openmoc.PERIODIC)
        p2y.setBoundaryType(openmoc.REFLECTIVE)
        p1z.setBoundaryType(openmoc.PERIODIC)
        p2z.setBoundaryType(openmoc.VACUUM)

        # Define halfspaces
        h1x = openmoc.Halfspace(-1, p1x)
        h2x = openmoc.Halfspace(+1, p2x)
        h1y = openmoc.Halfspace(-1, p1y)
        h2y = openmoc.Halfspace(+1, p2y)
        h1z = openmoc.Halfspace(-1, p1z)
        h2z = openmoc.Halfspace(+1, p2z)

        # Add halfspaces
        intersection.addNode(h1x)
        intersection.addNode(h2x)
        intersection.addNode(h1y)
        intersection.addNode(h2y)
        intersection.addNode(h1z)
        intersection.addNode(h2z)

        # Test getMaxXYZ + MaxXYZboundary
        self.assertEqual(intersection.getMaxX(), 3)
        self.assertEqual(intersection.getMaxXBoundaryType(), openmoc.VACUUM)
        self.assertEqual(intersection.getMaxY(), 1)
        self.assertEqual(intersection.getMaxYBoundaryType(), openmoc.PERIODIC)
        self.assertEqual(intersection.getMaxZ(), 8)
        self.assertEqual(intersection.getMaxZBoundaryType(), openmoc.PERIODIC)

        # Test getMinXYZ + MinXYZboundary
        self.assertEqual(intersection.getMinX(), -2)
        self.assertEqual(intersection.getMinXBoundaryType(), openmoc.REFLECTIVE)
        self.assertEqual(intersection.getMinY(), -0.5)
        self.assertEqual(intersection.getMinYBoundaryType(), openmoc.REFLECTIVE)
        self.assertEqual(intersection.getMinZ(), -4)
        self.assertEqual(intersection.getMinZBoundaryType(), openmoc.VACUUM)

    def test_complement(self):

        complement = openmoc.Complement()
        intersection = openmoc.Intersection()

        # Define surfaces
        p1x = openmoc.XPlane(x = 3)
        p2x = openmoc.XPlane(x = -2)
        p1y = openmoc.YPlane(y = 1)
        p2y = openmoc.YPlane(y = -0.5)
        p1z = openmoc.ZPlane(z = 8)
        p2z = openmoc.ZPlane(z = -4)

        # Define boundary types
        p1x.setBoundaryType(openmoc.VACUUM)
        p2x.setBoundaryType(openmoc.REFLECTIVE)
        p1y.setBoundaryType(openmoc.PERIODIC)
        p2y.setBoundaryType(openmoc.REFLECTIVE)
        p1z.setBoundaryType(openmoc.PERIODIC)
        p2z.setBoundaryType(openmoc.VACUUM)

        # Define halfspaces
        h1x = openmoc.Halfspace(-1, p1x)
        h2x = openmoc.Halfspace(+1, p2x)
        h1y = openmoc.Halfspace(-1, p1y)
        h2y = openmoc.Halfspace(+1, p2y)
        h1z = openmoc.Halfspace(-1, p1z)
        h2z = openmoc.Halfspace(+1, p2z)

        # Add halfspaces
        intersection.addNode(h1x)
        intersection.addNode(h2x)
        intersection.addNode(h1y)
        intersection.addNode(h2y)
        intersection.addNode(h1z)
        intersection.addNode(h2z)
        complement.addNode(intersection)

        # Test getMaxXYZ + MaxXYZboundary
        #FIXME Implement boundary type getter for complement regions
        self.assertEqual(complement.getMaxX(), 3)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMaxXBoundaryType(), openmoc.VACUUM)
        self.assertEqual(complement.getMaxY(), 1)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMaxYBoundaryType(), openmoc.PERIODIC)
        self.assertEqual(complement.getMaxZ(), 8)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMaxZBoundaryType(), openmoc.PERIODIC)

        # Test getMinXYZ + MinXYZboundary
        self.assertEqual(complement.getMinX(), -2)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMinXBoundaryType(), openmoc.REFLECTIVE)
        self.assertEqual(complement.getMinY(), -0.5)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMinYBoundaryType(), openmoc.REFLECTIVE)
        self.assertEqual(complement.getMinZ(), -4)
        with self.assertRaises(Exception): self.assertEqual(
             complement.getMinZBoundaryType(), openmoc.VACUUM)

if __name__ == '__main__':
    unittest.main()
