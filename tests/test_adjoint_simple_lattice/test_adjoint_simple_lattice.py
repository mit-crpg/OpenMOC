#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc


class SimpleLatticeAdjointTestHarness(TestHarness):
    """An adjoint eigenvalue calculation for a 4x4 lattice with 7-group C5G7
    cross section data."""

    def __init__(self):
        super(SimpleLatticeAdjointTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self.calculation_mode = openmoc.ADJOINT

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=True):
        """Digest info in the solver and return hash as a string."""
        return super(SimpleLatticeAdjointTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)


if __name__ == '__main__':
    harness = SimpleLatticeAdjointTestHarness()
    harness.main()
