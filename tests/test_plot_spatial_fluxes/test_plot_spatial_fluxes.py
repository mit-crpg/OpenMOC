#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import SimpleLatticeInput
from openmoc.plotter import plot_spatial_fluxes


class PlotSpatialFluxesTestHarness(PlottingTestHarness):
    """Test spatial flux plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotSpatialFluxesTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

    def _run_openmoc(self):
        """Run OpenMOC and plot the spatial fluxes in the geometry."""

        # Run an eigenvalue calculation
        super(PlotSpatialFluxesTestHarness, self)._run_openmoc()

        # Specify energy groups for which to plot the spatial flux
        energy_groups = [1, 3, 5, 7]

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100,
                       get_figure=True, energy_groups=energy_groups))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, get_figure=True,
                       xlim=(0., 2.), ylim=(0., 2.), 
                       energy_groups=energy_groups))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, get_figure=True,
                       energy_groups=energy_groups, library='pil'))


if __name__ == '__main__':
    harness = PlotSpatialFluxesTestHarness()
    harness.main()
