#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import GridInput
import openmoc


class MultiSimMaterialsCmfdTestHarness(TestHarness):
    """A multi-simulation eigenvalue calculation with different materials in
     a grid geometry with CMFD."""

    def __init__(self):
        super(MultiSimMaterialsCmfdTestHarness, self).__init__()
        self.input_set = GridInput()
        self.num_simulations = 3
        self.num_iters = []
        self.keffs = []

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(MultiSimMaterialsCmfdTestHarness, self)._create_geometry()

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
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

            # Exchange all of the materials for their clones
            for cell_id in cells:
                material = cells[cell_id].getFillMaterial()
                material = material.clone()
                cells[cell_id].setFill(material)

            # Run eigenvalue calculation and store the results
            super(MultiSimMaterialsCmfdTestHarness, self)._run_openmoc()
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
    harness = MultiSimMaterialsCmfdTestHarness()
    harness.main()
