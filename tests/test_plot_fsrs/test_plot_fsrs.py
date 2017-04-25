#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import SimpleLatticeInput
from openmoc.plotter import plot_flat_source_regions


class PlotFSRsTestHarness(PlottingTestHarness):
    """Test flat source region plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotFSRsTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

    def _run_openmoc(self):
        """Plot the flat source regions in the geometry."""

        # Run an eigenvalue calculation to setup FSR centroids
        super(PlotFSRsTestHarness, self)._run_openmoc()

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100, 
                       get_figure=True))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100, 
                       get_figure=True, xlim=(0., 2.), ylim=(0., 2.)))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100, 
                       get_figure=True, centroids=True, marker_size=3))
        self.figures.append(
            plot_flat_source_regions(self.input_set.geometry, gridsize=100, 
                       get_figure=True, centroids=True, library='pil'))
 

if __name__ == '__main__':
    harness = PlotFSRsTestHarness()
    harness.main()
