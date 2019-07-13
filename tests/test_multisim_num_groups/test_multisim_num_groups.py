#!/usr/bin/env python

import os
import sys

import numpy as np

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import HomInfMedInput
import openmoc


class MultiSimNumGroupsTestHarness(MultiSimTestHarness):
    """A multi-simulation eigenvalue calculation for a homogeneous infinite
    medium with 1-group and 2-group data."""

    def __init__(self):
        super(MultiSimNumGroupsTestHarness, self).__init__()
        self.input_set = HomInfMedInput()
        self.num_simulations = 1

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with 1- and 2-group
        cross sections in a homogeneous infinite medium."""

        # Extract infinite medium material from input set
        material = self.input_set.materials['infinite medium']

        # Setup 1-group multi-group cross sections
        material.setName('1-group infinite medium')
        material.setNumEnergyGroups(1)
        material.setNuSigmaF(np.array([0.0994076580]))
        material.setSigmaS(np.array([0.383259177]))
        material.setChi(np.array([1.0]))
        material.setSigmaT(np.array([0.452648699]))

        # Run eigenvalue calculation and store the results
        super(MultiSimNumGroupsTestHarness, self)._run_openmoc()

        # Setup 2-group multi-group cross sections
        material.setName('2-group infinite medium')
        material.setNumEnergyGroups(2)
        material.setNuSigmaF(np.array([0.0015, 0.325]))
        material.setSigmaS(np.array([0.1, 0.117, 0., 1.42]))
        material.setChi(np.array([1.0, 0.0]))
        material.setSigmaT(np.array([0.2208, 1.604]))

        # Run eigenvalue calculation and store the results
        super(MultiSimNumGroupsTestHarness, self)._run_openmoc()


if __name__ == '__main__':
    harness = MultiSimNumGroupsTestHarness()
    harness.main()
