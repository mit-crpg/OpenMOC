#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import PinCellInput
from openmoc.plotter import plot_tracks


class PlotTracksTestHarness(PlottingTestHarness):
    """Test track plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotTracksTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Plot the tracks."""

        # Run an eigenvalue calculation to setup track generator
        super(PlotTracksTestHarness, self)._run_openmoc()

        # Create Matplotlib Figures for the tracks
        self.figures = [plot_tracks(self.track_generator, get_figure=True),
             plot_tracks(self.track_generator, get_figure=True, plot_3D=True)]


if __name__ == '__main__':
    harness = PlotTracksTestHarness()
    harness.main()
