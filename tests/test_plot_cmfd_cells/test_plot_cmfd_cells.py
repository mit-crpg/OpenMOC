#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import PwrAssemblyInput
import openmoc
from openmoc.plotter import plot_cmfd_cells


class PlotCmfdCellsTestHarness(PlottingTestHarness):
    """Test CMFD cells plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotCmfdCellsTestHarness, self).__init__()
        self.cmfd = None
        self.input_set = PwrAssemblyInput()

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(PlotCmfdCellsTestHarness, self)._create_geometry()

        # Initialize CMFD
        self.cmfd = openmoc.Cmfd()
        self.cmfd.setSORRelaxationFactor(1.5)
        self.cmfd.setLatticeStructure(17,17)
        self.cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
        self.cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(self.cmfd)

    def _run_openmoc(self):
        """Plot the flat source regions in the geometry."""

        # Run an eigenvalue calculation to setup CMFD cells
        super(PlotCmfdCellsTestHarness, self)._run_openmoc()

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd,
                            gridsize=100, get_figure=True))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            get_figure=True, xlim=(0., 2.), ylim=(0., 2.)))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            get_figure=True, library='pil'))
 

if __name__ == '__main__':
    harness = PlotCmfdCellsTestHarness()
    harness.main()
