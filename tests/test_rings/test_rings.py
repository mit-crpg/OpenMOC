#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class RingTestHarness(TestHarness):
    """Tests cell radial discretization."""

    def __init__(self):
        super(RingTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def run_openmoc(self):
        """Discretize fuel and moderator into rings."""

        self.input_set.create_materials()
        self.input_set.create_geometry()

        # Create different rings for each Cell
        cells = self.input_set.geometry.getAllMaterialCells()
        for i, cell_id in enumerate(cells):
            cells[cell_id].setNumRings(i*2 + 3)

        self.input_set.geometry.initializeFlatSourceRegions()
        self.create_trackgenerator()
        self.generate_tracks()

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Digest info from geometry and return as a string."""

        return super(RingTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, num_fsrs=num_fsrs,
                num_segments=num_segments, num_tracks=num_tracks)


if __name__ == '__main__':
    harness = RingTestHarness()
    harness.main()
