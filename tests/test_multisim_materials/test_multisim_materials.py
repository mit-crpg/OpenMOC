#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput
import openmoc


class MultiSimMaterialsTestHarness(TestHarness):
    """A multi-simulation eigenvalue calculation with different materials in
    a pin cell problem with C5G7 7-group cross section data."""

    def __init__(self):
        super(MultiSimMaterialsTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self.num_simulations = 3
        self.num_iters = []
        self.keffs = []

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations with different
        materials."""

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
            super(MultiSimMaterialsTestHarness, self)._run_openmoc()
            self.num_iters.append(self.solver.getNumIterations())
            self.keffs.append(self.solver.getKeff())

    def _get_results(self, num_iterations=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Return eigenvalues from each simulation into a string."""

        # Write out the iteration count and eigenvalues from each simulation
        outstr = ''
        for num_iters, keff in zip(self.num_iters, self.keffs):
            outstr += 'Iters: {0}\tkeff: {1:12.5E}\n'.format(num_iters, keff)

        return outstr


if __name__ == '__main__':
    harness = MultiSimMaterialsTestHarness()
    harness.main()
