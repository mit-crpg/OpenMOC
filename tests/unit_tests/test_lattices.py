import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestLattice(unittest.TestCase):

    def setUp(self):
        self.lattice = openmoc.Lattice(id=12, name="test lattice")
        u1 = openmoc.Universe(name='Universe 1')
        self.universe = u1
        u2 = openmoc.Universe(name='Universe 2')
        u3 = openmoc.Universe(name='Universe 3')
        self.lattice.setUniverses([[[u1, u2, u1, u2],
                             [u2, u3, u2, u3],
                             [u1, u2, u1, u2],
                             [u2, u3, u2, u3]],
                            [[u1, u2, u1, u2],
                             [u2, u3, u2, u3],
                             [u1, u2, u1, u2],
                             [u2, u3, u2, u3]]])
        self.cell = openmoc.Cell(name="test cell")
        self.lattice.addCell(self.cell)

    def test_getters(self):

        self.assertEqual(self.lattice.getId(), 12)
        self.assertEqual(self.lattice.getName(), "test lattice")
        self.assertEqual(self.lattice.getUid(), 13)

    def test_ids(self):

        openmoc.reset_universe_id()
        lattice_2 = openmoc.Lattice()
        self.assertEqual(lattice_2.getId(), 1000016)
        openmoc.maximize_universe_id(10000000)
        lattice_3 = openmoc.Lattice()
        self.assertEqual(lattice_3.getId(), 10000000)

    def test_cells(self):

        self.assertEqual(self.lattice.getNumCells(), 1)
        self.assertEqual(self.lattice.getCell(self.cell.getId()).getId(),
                         self.cell.getId())
        self.lattice.removeCell(self.cell)
        self.assertEqual(self.lattice.getNumCells(), 0)
        self.lattice.printString()

    def test_clone(self):

        self.cell.setFill(openmoc.Material())
        clone = self.lattice.clone()
        self.assertEqual(self.lattice.getUniverse(0,0,0).getName(),
                         "Universe 2")

    def test_fissionable(self):

        self.assertEqual(self.lattice.isFissionable(), False)
        self.lattice.setFissionability(True)
        self.assertEqual(self.lattice.isFissionable(), True)

    def test_universes(self):

        u4 = openmoc.Universe(name="Universe 4")
        self.lattice.updateUniverse(3,3,1,u4)
        self.assertEqual(self.lattice.getUniverse(3,3,1).getName(),
                         "Universe 4")

        self.lattice.removeUniverse(self.universe)
        self.assertEqual(self.lattice.getUniverse(0,1,0), None)

if __name__ == '__main__':
    unittest.main()
