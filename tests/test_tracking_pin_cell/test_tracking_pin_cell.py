#!/usr/bin/env python

import os
import sys
import math
from collections import OrderedDict
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TrackingTestHarness
from input_set import PinCellInput

import openmoc

class TrackingPinCellTestHarness(TrackingTestHarness):
    """Tests tracking over a pin cell geometry."""

    def __init__(self):
        super(TrackingPinCellTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _setup(self):
        """Initialize the materials, geometry, and tracks."""
        super(TrackingPinCellTestHarness, self)._create_geometry()

        # Initialize track objects
        tracks = OrderedDict()
        tracks['Diagonal Track'] = openmoc.Track()
        tracks['Tangent Track'] = openmoc.Track()
        tracks['Nudged Tangent Track'] = openmoc.Track()
        tracks['Horizontal Track'] = openmoc.Track()
        tracks['Vertical Track'] = openmoc.Track()
        tracks['Reverse Diagonal Track'] = openmoc.Track()

        # Set track trajectories and locations
        tracks['Diagonal Track'].setValues(-2, -2, 0, 2, 2, 0, math.atan(1))
        offset = math.sqrt(2) - 2
        tracks['Tangent Track'].setValues(offset, -2, 0, 2, -offset,\
                                                  0, math.atan(1))
        offset -= 1e-6
        tracks['Nudged Tangent Track'].setValues(offset, -2, 0, 2, -offset,\
                                                  0, math.atan(1))
        tracks['Horizontal Track'].setValues(-2, 0, 0, 2, 0, 0, 0)
        tracks['Vertical Track'].setValues(0, -2, 0, 0, 2, 0, math.pi/2)
        tracks['Reverse Diagonal Track'].setValues(2, 2, 0, -2, -2, 0,\
                                                   math.pi + math.atan(1))

        # Set the tracks dictionary in the parent class
        self.tracks = tracks

if __name__ == '__main__':
    harness = TrackingPinCellTestHarness()
    harness.main()
