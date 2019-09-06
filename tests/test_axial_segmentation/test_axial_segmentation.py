#!/usr/bin/env python

import os
import sys
import numpy as np
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import AxialExtendedInput


class AxialSegmentationTestHarness(TestHarness):
    """An eigenvalue calculation in a 3D lattice with 7-group C5G7 data."""

    def __init__(self):
        super(AxialSegmentationTestHarness, self).__init__()
        self.input_set = AxialExtendedInput()
        self.num_polar = 2
        self.azim_spacing = 0.24
        self.z_spacing = 0.9


    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(openmoc.OTF_TRACKS)

        # Add segmentation zones, missing some intentionally
        seg_zones = np.array([0., 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20])
        self.track_generator.setSegmentationZones(seg_zones);


    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Need to use more than 1 thread, and lose reproducibility of FSR
        # numbering, in order to have temporary tracks and segments array of
        # the correct size for the multi-threaded solver.
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()


    def _run_openmoc(self):
        self.solver.computeEigenvalue(max_iters=30)

    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver"""
        return super(AxialSegmentationTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = AxialSegmentationTestHarness()
    harness.main()
