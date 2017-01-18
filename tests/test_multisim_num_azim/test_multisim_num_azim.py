#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import PinCellInput
import openmoc


class MultiSimNumAzimTestHarness(MultiSimTestHarness):
    """A multi-simulation eigenvalue calculation with varying azimuthal angle
    counts for a simple pin cell with 7-group C5G7 cross section data."""

    def __init__(self):
        super(MultiSimNumAzimTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self.num_simulations = 1
        self.num_tracks = []
        self.num_segments = []

    def _setup(self):
        """Build materials, geometry and dummy track generator."""
        self._create_geometry()
        self._create_trackgenerator()
        self._create_solver()

    def _create_solver(self):
        """Instantiate a CPUSolver."""
        self.solver = openmoc.CPUSolver()
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations."""

        for num_azim in [4, 8, 16]:

            # Generate tracks
            self.track_generator.setNumAzim(num_azim)
            super(MultiSimNumAzimTestHarness, self)._generate_tracks()

            # Assign TrackGenerator to Solver and run eigenvalue calculation
            self.solver.setTrackGenerator(self.track_generator)
            super(MultiSimNumAzimTestHarness, self)._run_openmoc()

            # Store results
            self.num_tracks.append(self.track_generator.getNumTracks())
            self.num_segments.append(self.track_generator.getNumSegments())

    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=True, num_segments=True,
                     hash_output=False):
        """Return track and segment counts and eigenvalues from each
        simulation as a string."""

        outstr = ''

        # Write out the numbers of tracks and segments
        for num_tracks, num_segments in zip(self.num_tracks, self.num_segments):
            outstr += '# tracks: {0}\t'.format(num_tracks)
            outstr += '# segments: {0}\n'.format(num_segments)

        # Write out the iteration count and eigenvalues from each simulation
        for num_iters, keff in zip(self.num_iters, self.keffs):
            outstr += 'Iters: {0}\tkeff: {1:12.5E}\n'.format(num_iters, keff)

        return outstr


if __name__ == '__main__':
    harness = MultiSimNumAzimTestHarness()
    harness.main()
