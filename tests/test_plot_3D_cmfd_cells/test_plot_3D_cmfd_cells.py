#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import NonUniformLatticeInput
import openmoc
from openmoc.plotter import plot_cmfd_cells


class PlotCmfdCellsTestHarness(PlottingTestHarness):
    """Test CMFD cells plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotCmfdCellsTestHarness, self).__init__()
        self.cmfd = None
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
        super(PlotCmfdCellsTestHarness, self)._create_solver()
        # Use only 1 thread for FSR numbering reproducibility
        # and for OTF ray tracing
        self.solver.setNumThreads(1)

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(PlotCmfdCellsTestHarness, self)._create_geometry()

        # Initialize CMFD
        self.cmfd = openmoc.Cmfd()
        self.cmfd.setSORRelaxationFactor(1.5)
        self.cmfd.setLatticeStructure(4,4,3)
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
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, offset=0.1,
                            plane='xy', gridsize=100, get_figure=True))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            offset=0.1, plane='yz', get_figure=True,
                            ylim=(0., 2.), zlim=(0., 2.)))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            offset=0.1, plane='yz', get_figure=True,
                            library='pil'))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            offset=0.1, plane='yz', get_figure=True))
        self.figures.append(
            plot_cmfd_cells(self.input_set.geometry, self.cmfd, gridsize=100,
                            offset=0.1, plane='xz', get_figure=True))


if __name__ == '__main__':
    harness = PlotCmfdCellsTestHarness()
    harness.main()
