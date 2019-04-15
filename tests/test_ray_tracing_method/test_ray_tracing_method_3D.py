#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import SimpleLatticeInput


class RayTracingMethodTestHarness(TestHarness):
    """3D lattice eigenvalue calculation to test agreement of ray-tracers."""

    def __init__(self):
        super(RayTracingMethodTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(num_dimensions=3)
        self.num_polar = 4
        self.azim_spacing = 0.4
        self.z_spacing = 1.2
        self.tolerance = 1E-3

        # To store results
        self.num_iters = []
        self.keff = []
        self.fluxes = []
        self.num_fsrs = []
        self.num_tracks = []
        self.num_segments = []

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(RayTracingMethodTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(1.0)
        cmfd.setSORRelaxationFactor(1.5)
        cmfd.setLatticeStructure(4,4,4)
        cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
        cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _create_trackgenerator(self, segmentation_method=openmoc.OTF_TRACKS):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(segmentation_method)

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

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations."""

        for segmentation_method in [openmoc.OTF_TRACKS, openmoc.OTF_STACKS,
                                    openmoc.EXPLICIT_3D]:

            # Create a new geometry to reset problem
            super(RayTracingMethodTestHarness, self)._create_geometry()

            # Create track generator with chosen ray tracing method
            self._create_trackgenerator(segmentation_method)
            super(RayTracingMethodTestHarness, self)._generate_tracks()

            # Assign TrackGenerator to Solver and run eigenvalue calculation
            self.solver.setTrackGenerator(self.track_generator)
            super(RayTracingMethodTestHarness, self)._run_openmoc()

            # Store results
            self.num_iters.append(self.solver.getNumIterations())
            self.keff.append(self.solver.getKeff())
            out_fluxes = self.solver.getFluxes(7 * self.input_set.geometry.getNumFSRs())
            self.fluxes.append(out_fluxes)
            self.num_fsrs.append(self.input_set.geometry.getNumFSRs())
            self.num_tracks.append(self.track_generator.getNumTracks())
            self.num_segments.append(self.track_generator.getNumSegments())

    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=False):

        msg_otf = "OTF_TRACKS and OTF_STACKS results don't match"
        assert self.num_iters[0] == self.num_iters[1], msg_otf
        assert round(self.keff[0], 6) == round(self.keff[1], 6), msg_otf
        assert self.num_fsrs[0] == self.num_fsrs[1], msg_otf
        assert self.num_tracks[0] == self.num_tracks[1], msg_otf
        assert self.num_segments[0] == self.num_segments[1], msg_otf
        # Fluxes cannot be compared easily since FSR numberings don't match

        msg_exp = "EXPLICIT_3D and OTF_STACKS results don't match"
        assert self.num_iters[2] == self.num_iters[1], msg_exp
        assert round(self.keff[2], 6) == round(self.keff[1], 6), msg_exp
        assert self.num_fsrs[2] == self.num_fsrs[1], msg_exp
        assert self.num_tracks[2] == self.num_tracks[1], msg_exp
        assert self.num_segments[2] == self.num_segments[1], msg_exp

        """Digest info in the solver"""
        return super(RayTracingMethodTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = RayTracingMethodTestHarness()
    harness.main()
