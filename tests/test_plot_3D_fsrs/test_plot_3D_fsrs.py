#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import NonUniformLatticeInput
from openmoc.plotter import plot_flat_source_regions
import openmoc

class PlotFSRsTestHarness(PlottingTestHarness):
    """Test flat source region plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotFSRsTestHarness, self).__init__()
        self.input_set = NonUniformLatticeInput()
        self.num_polar = 4
        self.azim_spacing = 0.5
        self.z_spacing = 2.0
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
        super(PlotFSRsTestHarness, self)._create_solver()
        # Use only 1 thread for FSR numbering reproducibility
        self.solver.setNumThreads(1)

    def _run_openmoc(self):
        """Plot the flat source regions in the geometry."""

        # Run an eigenvalue calculation to setup FSR centroids
        super(PlotFSRsTestHarness, self)._run_openmoc()

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100,
                       offset=0.1, get_figure=True))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100,
                       offset=0.1, plane='yz', get_figure=True, ylim=(0., 2.),
                       zlim=(0., 2.)))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100,
                       offset=0.1, plane='yz', get_figure=True, centroids=True,
                       marker_size=3))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100,
                       offset=0.1, plane='xz', get_figure=True, centroids=True,
                       library='pil'))


if __name__ == '__main__':
    harness = PlotFSRsTestHarness()
    harness.main()
