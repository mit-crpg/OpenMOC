import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class SurfacesTestCases:

    # Class for shared tests between all surface classes
    class TestSurfaces(unittest.TestCase):

        def test_getters(self):

            self.assertEqual(self.surface.getId(), 12)
            self.assertEqual(self.surface.getName(), "test surface")
            self.assertEqual(self.surface.getBoundaryType(), openmoc.VACUUM)

        def test_setters(self):

            self.surface.setName("new name")
            self.surface.setBoundaryType(openmoc.REFLECTIVE)
            self.assertEqual(self.surface.getName(), "new name")
            self.assertEqual(self.surface.getBoundaryType(), openmoc.REFLECTIVE)

class ZCylinderTests(SurfacesTestCases.TestSurfaces):

    def setUp(self):
        self.surface = openmoc.ZCylinder(x=4, y=1, radius=10,
                                        name="test surface", id=12)
        self.surface.setBoundaryType(openmoc.VACUUM)

    def test_surfaceType(self):
        self.assertEqual(self.surface.getSurfaceType(), openmoc.ZCYLINDER)

    def test_bounds(self):

        self.assertEqual(self.surface.getMinX(-1), -6)
        self.assertEqual(self.surface.getMaxX(-1), 14)
        self.assertEqual(self.surface.getMinY(-1), -9)
        self.assertEqual(self.surface.getMaxY(-1), 11)
        self.assertLess(self.surface.getMinZ(-1), -1e20)
        self.assertGreater(self.surface.getMaxZ(-1), 1e20)

    def test_distance(self):

        point = openmoc.Point()

        # minDistance calls intersection()
        # No intersection
        point.setX(0)
        point.setY(-10)
        self.assertGreater(self.surface.getMinDistance(point, 0, np.pi/4), 1e20)

        # Tangeant intersection
        point.setX(4)
        point.setY(1 - self.surface.getRadius() * np.sqrt(2))
        self.assertEqual(self.surface.getMinDistance(point, np.pi/4, np.pi/4),
                         self.surface.getRadius() / (np.sqrt(2) / 2))

        # Two intersections, since inside the circle
        point.setX(4)
        point.setY(1)
        self.assertEqual(self.surface.getMinDistance(point, 0, np.pi/2),
                         self.surface.getRadius())

        # Test the vertical (in y) case as well
        # No intersection
        point.setX(20)
        point.setY(0)
        self.assertGreater(self.surface.getMinDistance(point, np.pi/2, np.pi/4), 1e20)

        # Tangeant intersection
        point.setX(-6)
        point.setY(20)
        self.assertAlmostEqual(self.surface.getMinDistance(point, 3*np.pi/2, np.pi/4),
                               19 * np.sqrt(2))

        # Two intersections, since inside the circle
        point.setX(4)
        point.setY(1)
        self.assertTrue(self.surface.getMinDistance(point, np.pi/2, np.pi/2),
                        self.surface.getRadius())

    def test_onsurface(self):

        point = openmoc.Point()
        for x,y in [(-6, 1), (14, 1), (4, -9), (4, 11)]:
            point.setX(x)
            point.setY(y)
            self.assertTrue(self.surface.isPointOnSurface(point))
            self.assertEqual(self.surface.evaluate(point), 0)

if __name__ == '__main__':
    unittest.main()
