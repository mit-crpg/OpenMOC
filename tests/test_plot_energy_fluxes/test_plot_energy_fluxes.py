#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import PinCellInput
from openmoc.plotter import plot_energy_fluxes


class PlotEnergyFluxesTestHarness(PlottingTestHarness):
    """Test spatial flux plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotEnergyFluxesTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Run OpenMOC and plot the spatial fluxes in the geometry."""

        # Run an eigenvalue calculation
        super(PlotEnergyFluxesTestHarness, self)._run_openmoc()

        # Extract the FSR count from the geometry
        num_fsrs = self.input_set.geometry.getNumFSRs()
        fsrs = tuple(range(num_fsrs))

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.extend(plot_energy_fluxes(self.solver, fsrs,
                    get_figure=True))
        self.figures.extend(plot_energy_fluxes(self.solver, fsrs,
                    get_figure=True, loglog=True))


if __name__ == '__main__':
    harness = PlotEnergyFluxesTestHarness()
    harness.main()
