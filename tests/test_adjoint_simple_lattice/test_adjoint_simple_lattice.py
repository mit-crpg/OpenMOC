#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc


class SimpleLatticeAdjointTestHarness(TestHarness):
    """An adjoint eigenvalue calculation for a 4x4 lattice with 7-group C5G7 
    cross section data."""

    def __init__(self):
        super(SimpleLatticeAdjointTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self.calculation_mode = openmoc.ADJOINT


if __name__ == '__main__':
    harness = SimpleLatticeAdjointTestHarness()
    harness.main()
