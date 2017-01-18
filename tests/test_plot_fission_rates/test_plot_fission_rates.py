#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import SimpleLatticeInput
from openmoc.plotter import plot_fission_rates


class PlotFissionRatesTestHarness(PlottingTestHarness):
    """Test fission rate plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotFissionRatesTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

    def _run_openmoc(self):
        """Run OpenMOC and plot the fission rates in the geometry."""

        # Run an eigenvalue calculation
        super(PlotFissionRatesTestHarness, self)._run_openmoc()

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100,
                       get_figure=True))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                       xlim=(0., 2.), ylim=(0., 2.)))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                       norm=False, transparent_zeros=False))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                        library='pil'))


if __name__ == '__main__':
    harness = PlotFissionRatesTestHarness()
    harness.main()
