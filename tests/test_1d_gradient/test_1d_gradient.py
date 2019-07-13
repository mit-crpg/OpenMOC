#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import HomInfMedInput
import openmoc


class OneDGradientTestHarness(TestHarness):
    """An eigenvalue calculation in a cube with vacuum BCs in x
    and reflective BCs in y with 2-group cross section data."""

    def _create_geometry(self):
        """Put VACUUM boundary conditions on left and right boundaries."""

        self.input_set.create_materials()
        self.input_set.create_geometry()

        # Get the root Cell
        cells = self.input_set.geometry.getAllCells()
        for cell_id in cells:
            cell = cells[cell_id]
            if cell.getName() == 'root cell':
                root_cell = cell

        # Apply VACUUM BCs on the min/max XPlane surfaces
        surfaces = root_cell.getSurfaces()
        for surface_id in surfaces:
            surface = surfaces[surface_id]._surface
            if len(surface.getName()) > 0:
                if surface.getName()[0] == 'x':
                    surface.setBoundaryType(openmoc.VACUUM)

    def __init__(self):
        super(OneDGradientTestHarness, self).__init__()
        self.input_set = HomInfMedInput()


if __name__ == '__main__':
    harness = OneDGradientTestHarness()
    harness.main()
