#!/usr/bin/env python

import os
import sys
import math
from collections import OrderedDict
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TrackingTestHarness
from input_set import LatticeGridInput

import openmoc

class TrackingLatticeGridTestHarness(TrackingTestHarness):
    """Tests tracking over a lattice geometry."""

    def __init__(self):
        super(TrackingLatticeGridTestHarness, self).__init__()
        self.input_set = LatticeGridInput()

    def _setup(self):
        """Initialize the materials, geometry, and tracks."""
        super(TrackingLatticeGridTestHarness, self)._create_geometry()

        # Initialize track objects
        tracks = OrderedDict()
        tracks['Diagonal Track'] = openmoc.Track()
        tracks['Nudged Diagonal Track'] = openmoc.Track()
        tracks['Horizontal Track'] = openmoc.Track()
        tracks['Vertical Track'] = openmoc.Track()
        tracks['Reverse Diagonal Track'] = openmoc.Track()

        # Set track trajectories and locations
        tracks['Diagonal Track'].setValues(-3, -3, 0, 3, 3, 0, math.atan(1))
        nudge = 1e-5
        tracks['Nudged Diagonal Track'].setValues(-3+nudge, -3, 0, 3, 3-nudge,\
                                                  0, math.atan(1))
        tracks['Horizontal Track'].setValues(-3, 0, 0, 3, 0, 0, 0)
        tracks['Vertical Track'].setValues(0, -3, 0, 0, 3, 0, math.pi/2)
        tracks['Reverse Diagonal Track'].setValues(3, 3, 0, -3, -3, 0,\
                                                   math.pi + math.atan(1))

        # Set the tracks dictionary in the parent class
        self.tracks = tracks

if __name__ == '__main__':
    harness = TrackingLatticeGridTestHarness()
    harness.main()
