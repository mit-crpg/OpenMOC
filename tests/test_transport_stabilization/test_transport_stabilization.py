#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc

class TransportStabilizationTestHarness(TestHarness):
    """An eigenvalue calculation for a 4x4 lattice with 7-group C5G7
    cross section data modified to present large negative in-scatter"""

    def __init__(self):
        super(TransportStabilizationTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(TransportStabilizationTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(0.7)
        cmfd.setLatticeStructure(17,17)
        cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
        cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _create_solver(self):
        """Instantiate a CPULSSolver."""
        self.solver = openmoc.CPULSSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with different
        stabilization techniques."""

        # Set inscatter to be negative
        material = self.input_set.materials['Water']
        material.setSigmaSByGroup(-1, 4, 4)

        stabilization_methods = [openmoc.DIAGONAL, openmoc.YAMAMOTO,
                                 openmoc.GLOBAL]

        for stab in stabilization_methods:

            self.solver.stabilizeTransport(0.4, stab)

            # Print negative sources once
            self.solver.printAllNegativeSources(stab == openmoc.DIAGONAL)

            # Run eigenvalue calculation and store the results
            self.num_simulations = 1
            super(TransportStabilizationTestHarness, self)._run_openmoc()

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=True):

        # Compare the negative sources distribution to the reference
        if (os.system("cmp k_negative_sources_iter_0 negative_sources_reference") != 0):
            py_printf('ERROR', "Prolongation ratios are different")

        """Digest info in the solver and return hash as a string."""
        return super(TransportStabilizationTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

    def _cleanup(self):
        """Delete the update ratio files."""
        super(TransportStabilizationTestHarness, self)._cleanup()
        os.remove("k_negative_sources_iter_0")


if __name__ == '__main__':
    harness = TransportStabilizationTestHarness()
    harness.main()
