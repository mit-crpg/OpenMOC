#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import PinCellInput
from openmoc.plotter import plot_quadrature


class PlotQuadratureTestHarness(PlottingTestHarness):
    """Test quadrature plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotQuadratureTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Plot the polar quadrature."""

        # Run an eigenvalue calculation to setup polar quadrature
        super(PlotQuadratureTestHarness, self)._run_openmoc()

        # Create Matplotlib Figures for the polar quadrature
        self.figures = [plot_quadrature(self.solver, get_figure=True)]


if __name__ == '__main__':
    harness = PlotQuadratureTestHarness()
    harness.main()
