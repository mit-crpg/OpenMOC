#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import GridInput

import openmoc

class TrackingGridTestHarness(TestHarness):
    """Tests tracking over a grid geometry."""

    def __init__(self):
        super(TrackingGridTestHarness, self).__init__()
        self.input_set = GridInput()
        self._result = ''

    def _setup(self):
        """Initialize the materials and geometry in the InputSet."""
        super(TrackingGridTestHarness, self)._create_geometry()

    def _segment_track(self, track, geometry):
        """Segments a given track over a given geometry and records the
           resulting segment information to a string"""

        # segmentize a track in a geometry, recording the segments in a string
        geometry.segmentize(track)
        num_segments = track.getNumSegments()
        info = 'Number of segments = ' + str(num_segments) + '\n'
        for i in range(num_segments):
            info += 'Segment ' + str(i) + ': '
            segment = track.getSegment(i)
            info += 'length=' + str(round(segment._length, 7)) + ', '
            info += 'FSR ID=' + str(segment._region_id) + ', '
            info += 'CMFD FWD=' + str(segment._cmfd_surface_fwd) + ', '
            info += 'CMFD BWD=' + str(segment._cmfd_surface_bwd) + ', '
            info += 'Material Name=' + str(segment._material.getName()) + ', '
            info += 'Material ID=' + str(segment._material.getId()) + '\n'
        track.clearSegments()
        return info

    def _run_openmoc(self):
        """Creates tracks over the geometry and segments them, saving the
           results in the _result string"""

        # initialize track objects
        diag_track = openmoc.Track()
        nudge_diag_track = openmoc.Track()
        hor_track = openmoc.Track()
        ver_track = openmoc.Track()
        rev_diag_track = openmoc.Track()
        geometry = self.input_set.geometry

        # set track trajectories and locations
        diag_track.setValues(-3, -3, 0, 3, 3, 0, math.atan(1))
        nudge = 1e-5
        nudge_diag_track.setValues(-3+nudge, -3, 0, 3, 3-nudge, 0, math.atan(1))
        hor_track.setValues(-3, 0, 0, 3, 0, 0, 0)
        ver_track.setValues(0, -3, 0, 0, 3, 0, math.pi/2)
        rev_diag_track.setValues(3, 3, 0, -3, -3, 0, math.pi + math.atan(1))

        # segmentize over the geometry
        self._result += 'Diagonal track\n'
        self._result += self._segment_track(diag_track, geometry)
        self._result += 'Nudged Diagonal track\n'
        self._result += self._segment_track(nudge_diag_track, geometry)
        self._result += 'Horizontal track\n'
        self._result += self._segment_track(hor_track, geometry)
        self._result += 'Vertical track\n'
        self._result += self._segment_track(ver_track, geometry)
        self._result += 'Reverse Diagonal track\n'
        self._result += self._segment_track(rev_diag_track, geometry)

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""
        return self._result

if __name__ == '__main__':
    harness = TrackingGridTestHarness()
    harness.main()
