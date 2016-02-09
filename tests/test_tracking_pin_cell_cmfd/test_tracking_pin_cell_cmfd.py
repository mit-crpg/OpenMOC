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

class TrackingPinCellCMFDTestHarness(TrackingTestHarness):
    """Tests tracking over a pin cell geometry with an overlaid CMFD mesh."""

    def __init__(self):
        super(TrackingPinCellCMFDTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _setup(self):
        """Initialize the materials, geometry, and tracks."""
        super(TrackingPinCellCMFDTestHarness, self)._create_geometry()

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

    def _run_openmoc(self):
        """Segment tracks over the geometry and save the result to a string"""

        # Segmentize over the geometry with a fine and coarse cmfd mesh
        for m in [3, 51]:

            # Overlay simple CMFD mesh
            self._result += '{0} x {0} CMFD mesh\n'.format(m)
            geometry = self.input_set.geometry
            cmfd = openmoc.Cmfd()
            cmfd.setLatticeStructure(m, m)
            geometry.setCmfd(cmfd)
            geometry.initializeCmfd()

            # Track over the composite geometry
            super(TrackingPinCellCMFDTestHarness, self)._run_openmoc()

if __name__ == '__main__':
    harness = TrackingPinCellCMFDTestHarness()
    harness.main()
