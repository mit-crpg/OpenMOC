import unittest
import numpy

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc


class TestCmfdParameters(unittest.TestCase):

    def test_cmfd_parameters(self):
        ''' Check default CMFD parameters and getters.'''

        cmfd = openmoc.Cmfd()

        # Check flux update
        self.assertEqual(cmfd.isFluxUpdateOn(), True)
        cmfd.setFluxUpdateOn(False)
        self.assertEqual(cmfd.isFluxUpdateOn(), False)

        # Check sigmaT rebalance
        self.assertEqual(cmfd.isSigmaTRebalanceOn(), False)
        cmfd.rebalanceSigmaT(True)
        self.assertEqual(cmfd.isSigmaTRebalanceOn(), True)

        # Check centroid update
        self.assertEqual(cmfd.isCentroidUpdateOn(), False)
        cmfd.setCentroidUpdateOn(True)
        self.assertEqual(cmfd.isCentroidUpdateOn(), True)

    def test_cmfd_group_structure(self):

        cmfd = openmoc.Cmfd()

        # CMFD number of groups getter test
        cmfd.setGroupStructure([[1, 2, 3], [4, 5]])
        self.assertEqual(cmfd.getNumCmfdGroups(), 2)

        # Check exception raising for non-monotonic group structure
        with self.assertRaises(Exception): cmfd.setGroupStructure([[1, 2, 6],
             [4, 5]])

        # Check exception raising for non-continuous group structure
        with self.assertRaises(Exception): cmfd.setGroupStructure([[1, 2, 3],
             [6, 7]])

if __name__ == '__main__':
    unittest.main()
