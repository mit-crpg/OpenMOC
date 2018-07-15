#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TrackingTestHarness
from input_set import GridInput

import openmoc

class TrackingGridTestHarness(TrackingTestHarness):
    """Tests tracking over a grid geometry."""

    def __init__(self):
        super(TrackingGridTestHarness, self).__init__()
        self.input_set = GridInput()

    def _setup(self):
        """Initialize the materials, geometry, and tracks."""
        super(TrackingGridTestHarness, self)._create_geometry()
        super(TrackingGridTestHarness, self)._create_trackgenerator()

        # Initialize track objects
        self.tracks['Diagonal Track'] = openmoc.Track()
        self.tracks['Nudged Diagonal Track'] = openmoc.Track()
        self.tracks['Horizontal Track'] = openmoc.Track()
        self.tracks['Vertical Track'] = openmoc.Track()
        self.tracks['Reverse Diagonal Track'] = openmoc.Track()

        # Set track trajectories and locations
        self.tracks['Diagonal Track'].setValues(-3, -3, 3, 3, math.atan(1))
        nudge = 1e-5
        self.tracks['Nudged Diagonal Track'].setValues(-3+nudge, -3, 3,\
                                                       3-nudge, math.atan(1))
        self.tracks['Horizontal Track'].setValues(-3, 0, 3, 0, 0)
        self.tracks['Vertical Track'].setValues(0, -3, 0, 3, math.pi/2)
        self.tracks['Reverse Diagonal Track'].setValues(3, 3, -3, -3,\
                                                        math.pi + math.atan(1))

if __name__ == '__main__':
    harness = TrackingGridTestHarness()
    harness.main()
