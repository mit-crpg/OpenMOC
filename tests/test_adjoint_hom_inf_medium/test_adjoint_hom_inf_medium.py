#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import HomInfMedInput
import openmoc


class HomInfMedAdjointTestHarness(TestHarness):
    """An adjoint eigenvalue calculation in a reflected cube with 2-group
    cross section data."""

    def __init__(self):
        super(HomInfMedAdjointTestHarness, self).__init__()
        self.input_set = HomInfMedInput()
        self.calculation_mode = openmoc.ADJOINT


if __name__ == '__main__':
    harness = HomInfMedAdjointTestHarness()
    harness.main()
