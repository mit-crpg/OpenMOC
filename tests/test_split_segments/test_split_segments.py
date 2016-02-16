#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class SplitSegmentsTestHarness(TestHarness):
    """An eigenvalue calculation in a pin cell with 7-group C5G7 data."""

    def __init__(self):
        super(SplitSegmentsTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue or fixed source calculation."""

        # Set a small max optical pathlength so segments are split
        self.solver.setMaxOpticalLength(0.5)

        super(SplitSegmentsTestHarness, self)._run_openmoc()

    def _get_results(self):
        """Digest info in the results and return as a string."""
        return super(SplitSegmentsTestHarness, self).\
            _get_results(num_segments=True)


if __name__ == '__main__':
    harness = SplitSegmentsTestHarness()
    harness.main()
