#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput

import openmoc
import math

class TrackingPinCellCMFDTestHarness(TestHarness):
    """Tests cell radial discretization."""

    def __init__(self):
        super(TrackingPinCellCMFDTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self._result = ''


    def _setup(self):
        """Initialize the materials and geometry in the InputSet."""
        super(TrackingPinCellCMFDTestHarness, self)._create_geometry()


    def _segment_track(self, track, geometry):

        # segmentize a track in a geometry, recording the segments in a string
        geometry.segmentize(track)
        num_segments = track.getNumSegments()
        info = 'Number of segments = ' + str(num_segments) + '\n'
        for i in range(num_segments):
            info += 'Segment ' + str(i) + ': '
            segment = track.getSegment(i)
            info += 'length=' + str(round(segment._length, 8)) + ', '
            info += 'FSR ID=' + str(segment._region_id) + ', '
            info += 'CMFD FWD=' + str(segment._cmfd_surface_fwd) + ', '
            info += 'CMFD BWD=' + str(segment._cmfd_surface_bwd) + ', '
            info += 'Material Name=' + str(segment._material.getName()) + ', '
            info += 'Material ID=' + str(segment._material.getId()) + '\n'
        track.clearSegments()
        return info

    def _run_openmoc(self):

        # initialize track objects
        diag_track = openmoc.Track()
        tan_track = openmoc.Track()
        nudge_tan_track = openmoc.Track()
        hor_track = openmoc.Track()
        ver_track = openmoc.Track()
        rev_diag_track = openmoc.Track()

        # set track trajectories and locations
        diag_track.setValues(-2, -2, 0, 2, 2, 0, math.atan(1))
        offset = math.sqrt(2) - 2
        tan_track.setValues(offset, -2, 0, 2, -offset, 0, math.atan(1))
        offset -= 1e-6
        nudge_tan_track.setValues(offset, -2, 0, 2, -offset, 0, math.atan(1))
        hor_track.setValues(-2, 0, 0, 2, 0, 0, 0.00)
        ver_track.setValues(0, -2, 0, 0, 2, 0, math.pi/2)
        rev_diag_track.setValues(2, 2, 0, -2, -2, 0, math.pi + math.atan(1))

        # segmentize over the geometry with a fine and coarse cmfd mesh
        for m in [3, 51]:

            # overlay simple CMFD mesh
            self._result += 'Pin cell with an overlaid {0} x {0} CMFD mesh\n'.format(m)
            geometry = self.input_set.geometry
            cmfd = openmoc.Cmfd()
            cmfd.setLatticeStructure(m, m)
            geometry.setCmfd(cmfd)
            geometry.initializeCmfd()

            # segmentize tracks over the geometry
            self._result += 'Diagonal track...\n'
            self._result += self._segment_track(diag_track, geometry)
            self._result += 'Tangent track...\n'
            self._result += self._segment_track(tan_track, geometry)
            self._result += 'Nudged Tangent track...\n'
            self._result += self._segment_track(nudge_tan_track, geometry)
            self._result += 'Horizontal track...\n'
            self._result += self._segment_track(hor_track, geometry)
            self._result += 'Vertical track...\n'
            self._result += self._segment_track(ver_track, geometry)
            self._result += 'Reverse Diagonal track...\n'
            self._result += self._segment_track(rev_diag_track, geometry)


    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""

        return self._result

if __name__ == '__main__':
    harness = TrackingPinCellCMFDTestHarness()
    harness.main()
