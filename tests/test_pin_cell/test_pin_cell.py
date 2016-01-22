#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput


class PinCellTestHarness(TestHarness):
    """An eigenvalue calculation in a pin cell with 7-group C5G7 data."""

    def __init__(self):
        super(PinCellTestHarness, self).__init__()
        self.input_set = PinCellInput()


if __name__ == '__main__':
    harness = PinCellTestHarness()
    harness.main()
