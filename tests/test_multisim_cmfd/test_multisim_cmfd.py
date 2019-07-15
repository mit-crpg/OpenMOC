#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import PwrAssemblyInput
import openmoc


class MultiSimCmfdTestHarness(MultiSimTestHarness):
    """A multi-simulation eigenvalue calculation with CMFD for a 17x17 
    lattice with 7-group C5G7 cross section data."""

    def __init__(self):
        super(MultiSimCmfdTestHarness, self).__init__()
        self.input_set = PwrAssemblyInput()
        self.max_iters = 5

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(MultiSimCmfdTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(1.0)
        cmfd.setSORRelaxationFactor(1.5)
        cmfd.setLatticeStructure(17,17)
        cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
        cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)


if __name__ == '__main__':
    harness = MultiSimCmfdTestHarness()
    harness.main()
