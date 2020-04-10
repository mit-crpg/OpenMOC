import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestCells(unittest.TestCase):

    def setUp(self):
        self.cell = openmoc.Cell(id=12, name="test cell")

    def test_getters(self):

        self.assertEqual(self.cell.getId(), 12)
        self.assertEqual(self.cell.getName(), "test cell")
        self.assertEqual(self.cell.getUid(), 0)

    def test_ids(self):

        openmoc.reset_cell_id()
        cell_2 = openmoc.Cell()
        self.assertEqual(cell_2.getId(), 1000000)
        openmoc.maximize_cell_id(10000000)
        cell_3 = openmoc.Cell()
        self.assertEqual(cell_3.getId(), 10000000)
        self.cell.printString()

    def test_volume(self):

        self.cell.setVolume(100)
        self.assertEqual(self.cell.getVolume(), 100)
        self.cell.incrementVolume(10)
        self.assertEqual(self.cell.getVolume(), 110)

    def test_instances(self):

        self.assertEqual(self.cell.getNumInstances(), 0)
        self.cell.setNumInstances(99)
        self.assertEqual(self.cell.getNumInstances(), 99)
        self.cell.incrementNumInstances()
        self.assertEqual(self.cell.getNumInstances(), 100)

    def test_rotations(self):

        # Test setting a rotation angle
        rotation = np.array([6, 2, 1])
        with self.assertRaises(Exception): self.cell.setRotation(np.array([0, 2.]))
        with self.assertRaises(Exception): self.cell.setRotation(rotation, "fake")
        material = openmoc.Material()
        self.cell.setFill(material)
        with self.assertRaises(Exception): self.cell.setRotation(rotation)

        universe = openmoc.Universe()
        self.cell.setFill(universe)
        self.cell.setRotation(rotation, "radians")

        # Test retrieving a rotation angle
        with self.assertRaises(Exception): self.cell.retrieveRotation(3, "fake")
        degrees = self.cell.retrieveRotation(3, "degrees")
        radians = self.cell.retrieveRotation(3, "radians")
        np.testing.assert_array_almost_equal(degrees, rotation / np.pi * 180, 10)
        np.testing.assert_array_almost_equal(radians, rotation, 10)


    def test_translations(self):

        # Test setting a rotation angle
        translation = np.array([6, 2, 1])
        with self.assertRaises(Exception): self.cell.setTranslation(np.array([0, 2.]))
        material = openmoc.Material()
        self.cell.setFill(material)
        with self.assertRaises(Exception): self.cell.setTranslation(translation)

        universe = openmoc.Universe()
        self.cell.setFill(universe)
        self.cell.setTranslation(translation)

        # Test retrieving a rotation angle
        output = self.cell.retrieveTranslation(3)
        np.testing.assert_array_almost_equal(translation, output, 10)

if __name__ == '__main__':
    unittest.main()
