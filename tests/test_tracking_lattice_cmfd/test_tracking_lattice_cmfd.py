#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TrackingTestHarness
from input_set import LatticeGridInput

import openmoc

class TrackingLatticeGridCMFDTestHarness(TrackingTestHarness):
    """Tests tracking over a lattice geometry with an overlaid CMFD mesh."""

    def __init__(self):
        super(TrackingLatticeGridCMFDTestHarness, self).__init__()
        self.input_set = LatticeGridInput()

    def _setup(self):
        """Initialize the materials, geometry, and tracks"""
        super(TrackingLatticeGridCMFDTestHarness, self)._create_geometry()
        super(TrackingLatticeGridCMFDTestHarness, self)._create_trackgenerator()

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
            super(TrackingLatticeGridCMFDTestHarness, self)._run_openmoc()

if __name__ == '__main__':
    harness = TrackingLatticeGridCMFDTestHarness()
    harness.main()
