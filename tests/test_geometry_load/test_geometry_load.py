#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import SimpleLatticeInput


class SimpleLatticeTestHarness(TestHarness):
    """An eigenvalue calculation in a 3D lattice with 7-group C5G7 data."""

    def __init__(self):
        super(SimpleLatticeTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _create_trackgenerator(self):
        """Dump and load the geometry then instantiate a TrackGenerator."""

        # Geometry should be dumped before FSRs are initialized
        # Dump the geometry
        self.input_set.geometry.dumpToFile("geometry_file.geo")

        # Get rid of the geometry
        self.input_set.geometry = None

        # Reload the geometry
        self.input_set.geometry = openmoc.Geometry()
        self.input_set.geometry.loadFromFile("geometry_file.geo")

        # Instantiate a TrackGenerator
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()

        self.track_generator = \
            openmoc.TrackGenerator(geometry, self.num_azim,
                                     self.spacing)

    def _create_solver(self):
        """Instantiate a CPULSSolver, to add the centroid calculations"""
        self.solver = openmoc.CPULSSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)


    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=False):
        """Digest info in the solver"""
        return super(SimpleLatticeTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = SimpleLatticeTestHarness()
    harness.main()
