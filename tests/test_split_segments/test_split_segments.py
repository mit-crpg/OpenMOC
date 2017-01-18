#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class SplitSegmentsTestHarness(TestHarness):
    """Test segment splitting based on max optical path length."""

    def __init__(self):
        super(SplitSegmentsTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _run_openmoc(self):
        """Set a small max optical path length to ensure segments are split."""

        # Set a small max optical path length so segments are split
        self.solver.setMaxOpticalLength(0.5)

        super(SplitSegmentsTestHarness, self)._run_openmoc()

    def _get_results(self):
        """Digest info in the results and return as a string."""
        return super(SplitSegmentsTestHarness, self).\
            _get_results(num_segments=True, fluxes=False, keff=False)


if __name__ == '__main__':
    harness = SplitSegmentsTestHarness()
    harness.main()
