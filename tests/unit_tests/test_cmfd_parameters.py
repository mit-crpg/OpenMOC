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

if __name__ == '__main__':
    unittest.main()
