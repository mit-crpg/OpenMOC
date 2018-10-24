#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import SimpleLatticeInput


class GeometryDumpTestHarness(TestHarness):
    """Test dumping a geometry to file."""

    def __init__(self):
        super(GeometryDumpTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(num_dimensions=3)


    def _create_solver(self):
        pass


    def _create_trackgenerator(self):
        pass


    def _generate_tracks(self):
        pass


    def _run_openmoc(self):

       # Dump geometry
       self.input_set.geometry.dumpToFile("geometry_file.geo")

       # NOTE : dumping the geometry before and after FSRs are initialized
       # leads to different results, as the rings and sectors create new cells.


    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest geometry file and return as a string."""

        # Find the geometry file
        filename = "geometry_file.geo"

        # Read the file into a string
        with open(filename, 'rb') as myfile:
            outstr = str(myfile.read())

        return outstr

if __name__ == '__main__':
    harness = GeometryDumpTestHarness()
    harness.main()
