#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import AxialExtendedInput
from openmoc.plotter import plot_segments
import openmoc


class PlotSegmentsTestHarness(PlottingTestHarness):
    """Test segment plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotSegmentsTestHarness, self).__init__()
        self.input_set = AxialExtendedInput(small=True)
        self.num_polar = 2
        self.azim_spacing = 1.35
        self.z_spacing = 10.0
        self.max_iters = 1

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(openmoc.OTF_STACKS)

    def _create_solver(self):
        super(PlotSegmentsTestHarness, self)._create_solver()
        # Use only 1 thread for FSR numbering reproducibility
        # and for OTF ray tracing
        self.solver.setNumThreads(1)

    def _run_openmoc(self):
        """Plot the tracks."""

        # Run an eigenvalue calculation to setup track generator
        super(PlotSegmentsTestHarness, self)._run_openmoc()

        # Create Matplotlib Figures for the tracks
        self.figures = \
             [plot_segments(self.track_generator, plot_3D=True, get_figure=True)]


if __name__ == '__main__':
    harness = PlotSegmentsTestHarness()
    harness.main()
