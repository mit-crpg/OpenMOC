#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput
import openmoc


class PinCellAdjointTestHarness(TestHarness):
    """An adjoint eigenvalue calculation in a pin cell with 7-group C5G7 
    cross section data."""

    def __init__(self):
        super(PinCellAdjointTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self.calculation_mode = openmoc.ADJOINT


if __name__ == '__main__':
    harness = PinCellAdjointTestHarness()
    harness.main()
