#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import NonUniformLatticeInput
from openmoc.plotter import plot_fission_rates
import openmoc


class PlotFissionRatesTestHarness(PlottingTestHarness):
    """Test fission rate plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotFissionRatesTestHarness, self).__init__()
        self.input_set = NonUniformLatticeInput()
        self.num_polar = 4
        self.azim_spacing = 0.5
        self.z_spacing = 2.0
        self.max_iters = 10

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(openmoc.OTF_STACKS)

    def _create_solver(self):
        super(PlotFissionRatesTestHarness, self)._create_solver()
        # Use only 1 thread for FSR numbering reproducibility
        self.solver.setNumThreads(1)

    def _run_openmoc(self):
        """Run OpenMOC and plot the fission rates in the geometry."""

        # Run an eigenvalue calculation
        super(PlotFissionRatesTestHarness, self)._run_openmoc()

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, offset=0.1, plane='xy',
                       get_figure=True))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                       offset=0.6, plane='yz', ylim=(0., 2.), zlim=(0., 2.)))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                       offset=0.6, plane='xz', norm=False, transparent_zeros=False))
        self.figures.append(
            plot_fission_rates(self.solver, gridsize=100, get_figure=True,
                        offset=1.2, plane='xz'))


if __name__ == '__main__':
    harness = PlotFissionRatesTestHarness()
    harness.main()
