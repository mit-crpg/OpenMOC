#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PwrAssemblyInput
import openmoc


class MultiSimCmfdTestHarness(TestHarness):
    """A multi-simulation eigenvalue calculation for a 17x17 lattice with 
    7-group C5G7 cross section data."""

    def __init__(self):
        super(MultiSimCmfdTestHarness, self).__init__()
        self.input_set = PwrAssemblyInput()
        self.keffs = []
        self.num_simulations = 3
        self.max_iters = 5

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(MultiSimCmfdTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setSORRelaxationFactor(1.5)
        cmfd.setLatticeStructure(17,17)
        cmfd.setGroupStructure([1,4,8])
        cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with CMFD."""

        for i in range(self.num_simulations):
            super(MultiSimCmfdTestHarness, self)._run_openmoc()
            
            self.keffs.append(self.solver.getKeff())

    def _get_results(self, num_iters=False, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Return eigenvalues from each simulation into a string."""

        outstr = ''

        # Write out the eigenvalues from each simulation
        for keff in self.keffs:
            outstr += 'keff: {0:12.5E}\n'.format(keff)

        return outstr


if __name__ == '__main__':
    harness = MultiSimCmfdTestHarness()
    harness.main()
