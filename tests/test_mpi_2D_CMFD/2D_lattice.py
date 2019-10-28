#!/usr/bin/env python

import os
import sys
from mpi4py import MPI

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc



class PwrAssemblyTestHarness(TestHarness):
    """An eigenvalue calculation for a 4x4 lattice with 7-group C5G7
    cross section data."""

    def __init__(self):
        super(PwrAssemblyTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(num_dimensions=2)
        self.spacing = 0.12
        self.max_iters = 1

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(PwrAssemblyTestHarness, self)._create_geometry()
        self.input_set.geometry.setDomainDecomposition(3, 2, 1, MPI.COMM_WORLD)

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setLatticeStructure(6,4)
        cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return hash as a string."""
        return super(PwrAssemblyTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)


if __name__ == '__main__':
    harness = PwrAssemblyTestHarness()
    harness.main()
