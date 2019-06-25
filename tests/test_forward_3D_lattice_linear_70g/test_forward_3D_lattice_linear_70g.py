#!/usr/bin/env python

import os
import sys
import numpy as np
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import SimpleLatticeInput


class SimpleLatticeTestHarness(TestHarness):
    """An eigenvalue calculation in a 3D lattice with 70g cross section data."""

    def __init__(self):
        super(SimpleLatticeTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(num_dimensions=3)
        self.num_polar = 2
        self.azim_spacing = 0.6
        self.z_spacing = 2.8
        self.tolerance = 5e-3


    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(openmoc.OTF_TRACKS)


    def _create_solver(self):
        """Instantiate a CPULSSolver."""
        self.solver = openmoc.CPULSSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)


    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Need to use more than 1 thread, and lose reproducibility of FSR
        # numbering, in order to have temporary tracks and segments array of
        # the correct size for the multi-threaded solver.
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()


    def _run_openmoc(self):
        # Extract UO2 material from input set, set 70-g cross sections
        material = self.input_set.materials['UO2']
        material.setNumEnergyGroups(70)
        material.setNuSigmaF(np.linspace(0, 1, 70) * 7)
        material.setSigmaS(np.linspace(0, 1, 4900) / 1000)
        material.setChi(1/70. * np.ones(70))
        material.setSigmaT(np.linspace(2, 3, 70))

        # Extract water material from input set, set 70-g cross sections
        material = self.input_set.materials['Water']
        material.setNumEnergyGroups(70)
        material.setNuSigmaF(np.linspace(0, 0, 70))
        material.setSigmaS(np.linspace(1, 2, 4900) / 1000)
        material.setChi(0/70. * np.ones(70))
        material.setSigmaT(np.linspace(3, 4, 70))

        self.solver.setConvergenceThreshold(5e-3)
        self.solver.computeEigenvalue(res_type=openmoc.FISSION_SOURCE)


    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver"""
        return super(SimpleLatticeTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = SimpleLatticeTestHarness()
    harness.main()
