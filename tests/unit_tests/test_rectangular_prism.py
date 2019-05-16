import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestRectangularPrism(unittest.TestCase):

    def test_constructor(self):
        rectprism = openmoc.RectangularPrism(1, 2, .5, 1.5, 3, 2.5)

        self.assertEqual(rectprism.getMinX(), 0)
        self.assertEqual(rectprism.getMinY(), .5)
        self.assertEqual(rectprism.getMinZ(), 1)
        self.assertEqual(rectprism.getMaxX(), 1)
        self.assertEqual(rectprism.getMaxY(), 2.5)
        self.assertEqual(rectprism.getMaxZ(), 4)

    def test_boundary_type_setter(self):
        rectprism = openmoc.RectangularPrism(1, 2, .5, 1.5, 3, 2.5)

        rectprism.setBoundaryType(openmoc.VACUUM)
        min_btx = rectprism.getMinXBoundaryType()

        self.assertEqual(min_btx, openmoc.VACUUM)

if __name__ == '__main__':
    unittest.main()