#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput

import openmoc

class NumFissionableFSRsTestHarness(TestHarness):
    """Tests num fissionable FSRs created for a simple pin-cell lattice."""

    def __init__(self):
        super(NumFissionableFSRsTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self._result = ''

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Return the number of fissionable FSRs."""

        num_fissionable_FSRs = self.solver.getNumFissionableFSRs()
        outstr = '# fissionable FSRs: {0}'.format(num_fissionable_FSRs)

        return outstr


if __name__ == '__main__':
    harness = NumFissionableFSRsTestHarness()
    harness.main()
