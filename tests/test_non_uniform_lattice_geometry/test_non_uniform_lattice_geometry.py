#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import NonUniformLatticeInput


class NonUniformLatticeTestHarness(TestHarness):
    """An eigenvalue calculation in a 3D Non-uniform Lattice 
       with 7-group C5G7 data."""

    def __init__(self):
        super(NonUniformLatticeTestHarness, self).__init__()
        self.input_set = NonUniformLatticeInput()
        self.num_polar = 4
        self.azim_spacing = self.spacing
        self.z_spacing = 0.5
        self.tolerance = 1E-4


    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(NonUniformLatticeTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(0.7)
        cmfd.setSORRelaxationFactor(1.0)
        cmfd.setWidths([[0.05,1.26,1.26,0.05], [0.05,1.26,1.26,0.05], [1.,1.5]])
        cmfd.setGroupStructure([[1,2,3],[4,5,6,7]])

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        
        quad = openmoc.EqualWeightPolarQuad()
        quad.setNumPolarAngles(self.num_polar)
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setQuadrature(quad)
        self.track_generator.setSegmentFormation(openmoc.OTF_STACKS)
        self.track_generator.setSegmentationZones([0.0,1.0,2.5])


    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Need to use more than 1 thread, and lose reproducibility of FSR
        # numbering, in order to have temporary tracks and segments array of
        # the correct size for the multi-threaded solver.
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()


    def _create_solver(self):
        """Instantiate a CPULSSolver."""
        self.solver = openmoc.CPULSSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)


    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=True, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver"""
        return super(NonUniformLatticeTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = NonUniformLatticeTestHarness()
    harness.main()
