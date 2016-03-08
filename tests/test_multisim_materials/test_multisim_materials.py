#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import PinCellInput
import openmoc


class MultiSimMaterialsTestHarness(MultiSimTestHarness):
    """A multi-simulation eigenvalue calculation with different materials in
    a pin cell problem with C5G7 7-group cross section data."""

    def __init__(self):
        super(MultiSimMaterialsTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with different
        materials."""

        for i in range(self.num_simulations):

            # Extract all of the Material-filled Cells in the Geometry
            cells = self.input_set.geometry.getAllMaterialCells()
            materials = self.input_set.geometry.getAllMaterials()

            # Exchange all of the Materials for their clones
            for cell_id in cells:
                material = cells[cell_id].getFillMaterial()
                clone = material.clone()
                cells[cell_id].setFill(clone)

            # Turn on SWIG flag to register old Materials 
            # with Python garbage collector
            for material_id in materials:
                materials[material_id].thisown = 1

            # Run eigenvalue calculation and store the results
            self.num_simulations = 1
            super(MultiSimMaterialsTestHarness, self)._run_openmoc()


if __name__ == '__main__':
    harness = MultiSimMaterialsTestHarness()
    harness.main()
