#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import HomInfMedInput


class HomInfMedTestHarness(TestHarness):
    """An eigenvalue calculation in a homogeneous infinite medium with 2-group
    cross section data."""

    def __init__(self):
        super(HomInfMedTestHarness, self).__init__()
        self.input_set = HomInfMedInput()


if __name__ == '__main__':
    harness = HomInfMedTestHarness()
    harness.main()
