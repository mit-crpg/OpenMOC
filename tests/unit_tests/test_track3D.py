import unittest
import numpy as np

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestTrack3D(unittest.TestCase):

    def test_setters_getters(self):

        # Set track3D components
        track = openmoc.Track3D()
        track.setTheta(0.5)
        track.setPolarIndex(1)
        track.setZIndex(2)
        track.setLZIndex(3)
        track.setCycleIndex(4)
        track.setCycleTrackIndex(5)
        track.setTrainIndex(6)

        # Use getter to check that setter worked
        self.assertEqual(track.getTheta(), 0.5)
        self.assertEqual(track.getPolarIndex(), 1)
        self.assertEqual(track.getZIndex(), 2)
        self.assertEqual(track.getLZIndex(), 3)
        self.assertEqual(track.getCycleIndex(), 4)
        self.assertEqual(track.getCycleTrackIndex(), 5)
        self.assertEqual(track.getTrainIndex(), 6)

if __name__ == '__main__':
    unittest.main()
