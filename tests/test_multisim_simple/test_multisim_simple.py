#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PwrAssemblyInput
import openmoc


class MultiSimSimpleTestHarness(TestHarness):
    """A multi-simulation eigenvalue calculation for a 17x17 lattice with 
    7-group C5G7 cross section data."""

    def __init__(self):
        super(MultiSimSimpleTestHarness, self).__init__()
        self.input_set = PwrAssemblyInput()
        self.keffs = []
        self.num_simulations = 3
        self.max_iters = 10

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations."""

        for i in range(self.num_simulations):
            super(MultiSimSimpleTestHarness, self)._run_openmoc()            
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
    harness = MultiSimSimpleTestHarness()
    harness.main()
