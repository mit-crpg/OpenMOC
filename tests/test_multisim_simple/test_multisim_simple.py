#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput
import openmoc


class MultiSimSimpleTestHarness(TestHarness):
    """A multi-simulation eigenvalue calculation for a 17x17 lattice with 
    7-group C5G7 cross section data."""

    def __init__(self):
        super(MultiSimSimpleTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self.num_simulations = 3
        self.num_iters = []
        self.keffs = []

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations."""

        for i in range(self.num_simulations):
            super(MultiSimSimpleTestHarness, self)._run_openmoc()            
            self.num_iters.append(self.solver.getNumIterations())
            self.keffs.append(self.solver.getKeff())

    def _get_results(self, num_iterations=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Return eigenvalues from each simulation into a string."""

        # Write out the iteration count and eigenvalues from each simulation
        outstr = ''
        for num_iters, keff in zip(self.num_iters, self.keffs):
            outstr += 'Iters: {0}\tkeff: {1:12.5E}\n'.format(num_iters, keff)

        return outstr


if __name__ == '__main__':
    harness = MultiSimSimpleTestHarness()
    harness.main()
