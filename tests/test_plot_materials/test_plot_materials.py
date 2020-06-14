#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import SimpleLatticeInput
from openmoc.plotter import plot_materials


class PlotMaterialsTestHarness(PlottingTestHarness):
    """Test material plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotMaterialsTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self.input_set.dimensions = 3

    def _run_openmoc(self):
        """Plot the materials in the geometry."""

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100,
                           get_figure=True))
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100,
                           offset=2.5, get_figure=True))
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100,
                           get_figure=True, xlim=(0., 2.), ylim=(0., 2.)))
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100,
                           get_figure=True, library='pil'))
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100, plane="yz",
                       get_figure=True))
        self.figures.append(
            plot_materials(self.input_set.geometry, gridsize=100, plane="xz",
                       get_figure=True))


if __name__ == '__main__':
    harness = PlotMaterialsTestHarness()
    harness.main()
