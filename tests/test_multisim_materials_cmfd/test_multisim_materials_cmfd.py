#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import GridInput
import openmoc


class MultiSimMaterialsCmfdTestHarness(MultiSimTestHarness):
    """A multi-simulation eigenvalue calculation with different materials in
     a grid geometry with CMFD."""

    def __init__(self):
        super(MultiSimMaterialsCmfdTestHarness, self).__init__()
        self.input_set = GridInput()

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(MultiSimMaterialsCmfdTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(1.0)
        cmfd.setSORRelaxationFactor(1.5)
        cmfd.setLatticeStructure(3,3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with CMFD with
        different materials."""

        for i in range(self.num_simulations):

            # Extract all of the material-filled cells in the geometry
            cells = self.input_set.geometry.getAllMaterialCells()
            materials = self.input_set.geometry.getAllMaterials()

            # Exchange all of the materials for their clones
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
            super(MultiSimMaterialsCmfdTestHarness, self)._run_openmoc()
            self.solver.printTimerReport()


if __name__ == '__main__':
    harness = MultiSimMaterialsCmfdTestHarness()
    harness.main()
