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

        # One intersections, from inside the cylinder
        point.setX(4)
        point.setY(1)
        self.assertEqual(self.surface.getMinDistance(point, 0, np.pi/2),
                         self.surface.getRadius())

        # Two intersections
        point.setX(-15)
        self.assertEqual(self.surface.getMinDistance(point, 0, np.pi/2), 9)

        # Test the vertical (in y) case as well, and test using coords
        # No intersection
        coords = openmoc.LocalCoords()
        coords.setX(20)
        coords.setY(0)
        coords.setPhi(np.pi/2)
        coords.setPolar(np.pi/4)
        self.assertGreater(self.surface.getMinDistance(coords), 1e20)

        # Tangeant intersection
        coords.setX(-6)
        coords.setY(20)
        coords.setPhi(3*np.pi/2)
        coords.setPolar(np.pi/4)
        self.assertAlmostEqual(self.surface.getMinDistance(coords),
                               19 * np.sqrt(2))

        # One intersection, from inside the circle
        coords.setX(4)
        coords.setY(1)
        coords.setPhi(np.pi/2)
        coords.setPolar(np.pi/2)
        self.assertTrue(self.surface.getMinDistance(coords),
                        self.surface.getRadius())

        # Two intersections
        coords.setY(15)
        coords.setPhi(np.pi/2)
        coords.setPolar(np.pi/2)
        self.assertTrue(self.surface.getMinDistance(coords), 6)

    def test_onsurface(self):

        point = openmoc.Point()
        for x,y in [(-6, 1), (14, 1), (4, -9), (4, 11)]:
            point.setX(x)
            point.setY(y)
            self.assertTrue(self.surface.isPointOnSurface(point))
            self.assertEqual(self.surface.evaluate(point), 0)


class PlaneTests(SurfacesTestCases.TestSurfaces):

    def setUp(self):
        self.surface = openmoc.Plane(A=2, B=1, C=2, D=0, name="test surface",
                                     id=12)
        self.surface.setBoundaryType(openmoc.VACUUM)

    def test_onsurface(self):

        point = openmoc.Point()
        point.setX(-8)
        point.setY(14)
        point.setZ(1)
        self.assertTrue(self.surface.isPointOnSurface(point))
        self.assertEqual(self.surface.evaluate(point), 0)


class XPlaneTests(PlaneTests):

    def setUp(self):
        self.surface = openmoc.XPlane(x=-8, name="test surface", id=12)
        self.surface.setBoundaryType(openmoc.VACUUM)

    def test_surfaceType(self):
        self.assertEqual(self.surface.getSurfaceType(), openmoc.XPLANE)

    def test_bounds(self):

        self.assertLess(self.surface.getMinZ(-1), -1e20)
        self.assertGreater(self.surface.getMaxZ(-1), 1e20)
        self.assertLess(self.surface.getMinY(-1), -1e20)
        self.assertGreater(self.surface.getMaxY(-1), 1e20)
        self.assertEqual(self.surface.getMinX(1), -8)
        self.assertEqual(self.surface.getMaxX(-1), -8)
        self.assertLess(self.surface.getMinX(-1), -1e20)
        self.assertGreater(self.surface.getMaxX(1), 1e20)

    def test_distance(self):

        point = openmoc.Point()

        # minDistance calls intersection()
        # No intersection
        point.setX(-10)
        point.setY(10)
        self.assertGreater(self.surface.getMinDistance(point, np.pi, 3*np.pi/4), 1e20)

        # Tangeant intersection : returns 0 intersection, as track parallel to the plane
        point.setX(-8)
        self.assertGreater(self.surface.getMinDistance(point, 0, np.pi/2), 1e20)

        # Intersection
        point.setX(-9)
        self.assertEqual(self.surface.getMinDistance(point, np.pi/4, np.pi/4), 2)


class YPlaneTests(PlaneTests):

    def setUp(self):
        self.surface = openmoc.YPlane(y=14, name="test surface", id=12)
        self.surface.setBoundaryType(openmoc.VACUUM)

    def test_surfaceType(self):
        self.assertEqual(self.surface.getSurfaceType(), openmoc.YPLANE)

    def test_bounds(self):

        self.assertLess(self.surface.getMinX(-1), -1e20)
        self.assertGreater(self.surface.getMaxX(-1), 1e20)
        self.assertLess(self.surface.getMinZ(-1), -1e20)
        self.assertGreater(self.surface.getMaxZ(-1), 1e20)
        self.assertEqual(self.surface.getMinY(1), 14)
        self.assertEqual(self.surface.getMaxY(-1), 14)
        self.assertLess(self.surface.getMinY(-1), -1e20)
        self.assertGreater(self.surface.getMaxY(1), 1e20)

    def test_distance(self):

        point = openmoc.Point()

        # minDistance calls intersection()
        # No intersection
        point.setX(-10)
        point.setY(10)
        self.assertGreater(self.surface.getMinDistance(point, np.pi, 3*np.pi/4), 1e20)

        # Tangeant intersection : returns 0 intersection, as track parallel to the plane
        point.setY(14)
        self.assertGreater(self.surface.getMinDistance(point, 0, np.pi/2), 1e20)

        # Intersection
        point.setY(13)
        self.assertEqual(self.surface.getMinDistance(point, np.pi/4, np.pi/4), 2)


class ZPlaneTests(SurfacesTestCases.TestSurfaces):

    def setUp(self):
        self.surface = openmoc.ZPlane(z=1, name="test surface", id=12)
        self.surface.setBoundaryType(openmoc.VACUUM)

    def test_surfaceType(self):
        self.assertEqual(self.surface.getSurfaceType(), openmoc.ZPLANE)

    def test_bounds(self):

        self.assertLess(self.surface.getMinX(-1), -1e20)
        self.assertGreater(self.surface.getMaxX(-1), 1e20)
        self.assertLess(self.surface.getMinY(-1), -1e20)
        self.assertGreater(self.surface.getMaxY(-1), 1e20)
        self.assertEqual(self.surface.getMinZ(1), 1)
        self.assertEqual(self.surface.getMaxZ(-1), 1)
        self.assertLess(self.surface.getMinZ(-1), -1e20)
        self.assertGreater(self.surface.getMaxZ(1), 1e20)

    def test_distance(self):

        point = openmoc.Point()

        # minDistance calls intersection()
        # No intersection, z of point is 0
        point.setX(4)
        point.setY(10)
        self.assertGreater(self.surface.getMinDistance(point, 0, 3*np.pi/4), 1e20)

        # Tangeant intersection : returns 0 intersection, as track parallel to the plane
        point.setZ(1)
        self.assertGreater(self.surface.getMinDistance(point, 0, np.pi/2), 1e20)

        # Intersection
        point.setZ(0)
        self.assertEqual(self.surface.getMinDistance(point, 0, np.pi/4), np.sqrt(2))

if __name__ == '__main__':
    unittest.main()
