#!/usr/bin/env python

import os
import sys
import numpy as np
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput
import openmoc

class PinCellTestHarness(TestHarness):
    """An eigenvalue calculation in a pin cell with 70-group cross sections."""

    def __init__(self):
        super(PinCellTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        # Extract UO2 material from input set, set 70-g cross sections
        material = self.input_set.materials['UO2']
        material.setNumEnergyGroups(70)
        material.setNuSigmaF(np.linspace(0, 1, 70) * 7)
        material.setSigmaS(np.linspace(0, 1, 4900) / 1000)
        material.setChi(1/70. * np.ones(70))
        material.setSigmaT(np.linspace(2, 3, 70))

        # Extract water material from input set, set 70-g cross sections
        material = self.input_set.materials['Water']
        material.setNumEnergyGroups(70)
        material.setNuSigmaF(np.linspace(0, 0, 70))
        material.setSigmaS(np.linspace(1, 2, 4900) / 1000)
        material.setChi(0/70. * np.ones(70))
        material.setSigmaT(np.linspace(3, 4, 70))

        self.solver.computeEigenvalue(res_type=openmoc.SCALAR_FLUX)

if __name__ == '__main__':
    harness = PinCellTestHarness()
    harness.main()
