#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import NonUniformLatticeInput
from openmoc.plotter import plot_spatial_fluxes
import openmoc


class PlotSpatialFluxesTestHarness(PlottingTestHarness):
    """Test spatial flux plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotSpatialFluxesTestHarness, self).__init__()
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
        super(PlotSpatialFluxesTestHarness, self)._create_solver()
        # Use only 1 thread for FSR numbering reproducibility
        self.solver.setNumThreads(1)

    def _run_openmoc(self):
        """Run OpenMOC and plot the spatial fluxes in the geometry."""

        # Run an eigenvalue calculation
        super(PlotSpatialFluxesTestHarness, self)._run_openmoc()

        # Specify energy groups for which to plot the spatial flux
        energy_groups = [1, 3, 5, 7]

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, offset=0.1,
                       get_figure=True, energy_groups=energy_groups))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, get_figure=True,
                       xlim=(0., 2.), ylim=(0., 2.), offset=0.1,
                       energy_groups=energy_groups))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, get_figure=True,
                       energy_groups=energy_groups, offset=0.1,
                       library='pil'))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, offset=0.1,
                       plane='yz', get_figure=True, energy_groups=energy_groups))
        self.figures.extend(
            plot_spatial_fluxes(self.solver, gridsize=100, get_figure=True,
                       xlim=(0., 2.), ylim=(0., 2.), offset=0.1, plane='xz',
                       energy_groups=energy_groups))

if __name__ == '__main__':
    harness = PlotSpatialFluxesTestHarness()
    harness.main()
