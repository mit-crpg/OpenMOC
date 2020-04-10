import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestUniverse(unittest.TestCase):

    def setUp(self):
        self.universe = openmoc.Universe(id=12, name="test universe")
        self.cell = openmoc.Cell(name="test cell")
        self.universe.addCell(self.cell)

    def test_getters(self):

        self.assertEqual(self.universe.getId(), 12)
        self.assertEqual(self.universe.getName(), "test universe")
        self.assertEqual(self.universe.getUid(), 4)

    def test_ids(self):

        openmoc.reset_universe_id()
        universe_2 = openmoc.Universe()
        self.assertEqual(universe_2.getId(), 1000001)
        openmoc.maximize_universe_id(10000000)
        universe_3 = openmoc.Universe()
        self.assertEqual(universe_3.getId(), 10000000)

    def test_cells(self):

        self.assertEqual(self.universe.getNumCells(), 1)
        self.assertEqual(self.universe.getCell(self.cell.getId()).getId(),
                         self.cell.getId())
        self.universe.removeCell(self.cell)
        self.assertEqual(self.universe.getNumCells(), 0)
        self.universe.printString()

    def test_clone(self):

        self.cell.setFill(openmoc.Material())
        clone = self.universe.clone()
        self.assertEqual(self.universe.getCell(self.cell.getId()).getName(),
                         "test cell")

    def test_fissionable(self):

        self.assertEqual(self.universe.isFissionable(), False)
        self.universe.setFissionability(True)
        self.assertEqual(self.universe.isFissionable(), True)

if __name__ == '__main__':
    unittest.main()
