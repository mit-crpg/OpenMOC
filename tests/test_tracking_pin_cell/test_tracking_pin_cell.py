#!/usr/bin/env python

import os
import sys
import math
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
        super(TrackingPinCellTestHarness, self)._create_trackgenerator()

        # Initialize track objects
        self.tracks['Diagonal Track'] = openmoc.Track()
        self.tracks['Tangent Track'] = openmoc.Track()
        self.tracks['Nudged Tangent Track'] = openmoc.Track()
        self.tracks['Horizontal Track'] = openmoc.Track()
        self.tracks['Vertical Track'] = openmoc.Track()
        self.tracks['Reverse Diagonal Track'] = openmoc.Track()

        # Set track trajectories and locations
        self.tracks['Diagonal Track'].setValues(-2, -2, 2, 2, math.atan(1))
        offset = math.sqrt(2) - 2
        self.tracks['Tangent Track'].setValues(offset, -2, 2, -offset,\
                                               math.atan(1))
        offset -= 1e-6
        self.tracks['Nudged Tangent Track'].setValues(offset, -2, 2,\
                                                      -offset, math.atan(1))
        self.tracks['Horizontal Track'].setValues(-2, 0, 2, 0, 0)
        self.tracks['Vertical Track'].setValues(0, -2, 0, 2, math.pi/2)
        self.tracks['Reverse Diagonal Track'].setValues(2, 2, -2, -2,\
                                                        math.pi + math.atan(1))

if __name__ == '__main__':
    harness = TrackingPinCellTestHarness()
    harness.main()
