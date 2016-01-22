#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class SectorTestHarness(TestHarness):
    """Tests cell sectorization"""

    def __init__(self):
        super(SectorTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def run_openmoc(self):
        """Run an OpenMOC eigenvalue or fixed source calculation."""

        self.input_set.create_materials()
        self.input_set.create_geometry()

        # Create different sectors for each Cell
        cells = self.input_set.geometry.getAllMaterialCells()
        for i, cell_id in enumerate(cells):
            cells[cell_id].setNumSectors(i*2 + 3)

        self.create_trackgenerator()
        self.generate_tracks()

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, hash_output=False):
        """Digest info from geometry and return as a string."""

        return super(SectorTestHarness, self)._get_results(num_iters=num_iters,
                                                           keff=keff,
                                                           fluxes=fluxes,
                                                           num_fsrs=num_fsrs,
                                                           hash_output=hash_output)


if __name__ == '__main__':
    harness = SectorTestHarness()
    harness.main()
