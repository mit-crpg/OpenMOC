#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput

import openmoc

class TrackingPinCellTestHarness(TestHarness):
    """Tests tracking over a pin cell geometry."""

    def __init__(self):
        super(TrackingPinCellTestHarness, self).__init__()
        self.input_set = PinCellInput()
        self._result = ''

    def _setup(self):
        """Initialize the materials and geometry in the InputSet."""
        super(TrackingPinCellTestHarness, self)._create_geometry()

    def _segment_track(self, track, geometry):
        """Segments a given track over a given geometry and records the
           resulting segment information to a string"""

        # Segmentize a track in a geometry, recording the segments in a string
        geometry.segmentize(track)
        num_segments = track.getNumSegments()
        info = ' ' + str(num_segments) + '\n'
        for i in range(num_segments):
            info += str(i) + ': '
            segment = track.getSegment(i)
            info += str(round(segment._length, 8)) + ', '
            info += str(segment._region_id) + ', '
            info += str(segment._cmfd_surface_fwd) + ', '
            info += str(segment._cmfd_surface_bwd) + ', '
            info += str(segment._material.getName()) + ', '
            info += str(segment._material.getId()) + '\n'
        track.clearSegments()
        return info

    def _run_openmoc(self):
        """Creates tracks over the geometry and segments them, saving the
           results in the _result string"""

        # Initialize track objects
        diag_track = openmoc.Track()
        tan_track = openmoc.Track()
        nudge_tan_track = openmoc.Track()
        hor_track = openmoc.Track()
        ver_track = openmoc.Track()
        rev_diag_track = openmoc.Track()
        geometry = self.input_set.geometry

        # Set track trajectories and locations
        diag_track.setValues(-2, -2, 0, 2, 2, 0, math.atan(1))
        offset = math.sqrt(2) - 2
        tan_track.setValues(offset, -2, 0, 2, -offset, 0, math.atan(1))
        offset -= 1e-6
        nudge_tan_track.setValues(offset, -2, 0, 2, -offset, 0, math.atan(1))
        hor_track.setValues(-2, 0, 0, 2, 0, 0, 0.00)
        ver_track.setValues(0, -2, 0, 0, 2, 0, math.pi/2)
        rev_diag_track.setValues(2, 2, 0, -2, -2, 0, math.pi + math.atan(1))

        # Segmentize over the geometry
        self._result += 'Diagonal track'
        self._result += self._segment_track(diag_track, geometry)
        self._result += 'Tangent track'
        self._result += self._segment_track(tan_track, geometry)
        self._result += 'Nudged Tangent track'
        self._result += self._segment_track(nudge_tan_track, geometry)
        self._result += 'Horizontal track'
        self._result += self._segment_track(hor_track, geometry)
        self._result += 'Vertical track'
        self._result += self._segment_track(ver_track, geometry)
        self._result += 'Reverse Diagonal track'
        self._result += self._segment_track(rev_diag_track, geometry)

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""
        return self._result

if __name__ == '__main__':
    harness = TrackingPinCellTestHarness()
    harness.main()
