#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
import openmoc

try:
    import openmc.openmoc_compatible
    import openmc.mgxs
except:
    print("OpenMC could not be imported, it's required for loading MGXS files")


class LoadMGXSTestHarness(MultiSimTestHarness):
    """Load a variety of OpenMC MGXS libraries."""

    def __init__(self):
        super(LoadMGXSTestHarness, self).__init__()
        self.azim_spacing = 0.12
        self.max_iters = 10
        self.keffs = []

        self.mgxs_files = ['mgxs_isotropic',
                          'mgxs_transport_corrected',
                          'mgxs_consistent',
                          'mgxs_consistent_nuscatter',
                          'mgxs_materials',
                          'mgxs_angular_legendre']
        # mgxs_angular_histogram currently not supported
        # mgxs_nuclide should be redone with the latest version of openmc
        # mgxs by distribcell, universe and mesh also not supported

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _create_solver(self):
        pass

    def _run_openmoc(self):
        """Run an OpenMOC calculation with each library."""

        for mgxs_file in self.mgxs_files:

            # Initialize OpenMC multi-group cross section library for a pin cell
            mgxs_lib = openmc.mgxs.Library.load_from_file(filename=mgxs_file,
                                                          directory='.')

            # Create an OpenMOC Geometry from the OpenMOC Geometry
            openmoc_geometry = \
                openmc.openmoc_compatible.get_openmoc_geometry(mgxs_lib.geometry)

            # Load cross section data
            openmoc_materials = \
                openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

            # Initialize FSRs
            openmoc_geometry.initializeFlatSourceRegions()

            # Initialize an OpenMOC TrackGenerator
            track_generator = openmoc.TrackGenerator(
                openmoc_geometry, self.num_azim, self.azim_spacing)
            track_generator.generateTracks()

            # Initialize an OpenMOC Solver
            self.solver = openmoc.CPUSolver(track_generator)
            self.solver.setConvergenceThreshold(self.tolerance)
            self.solver.setNumThreads(self.num_threads)

            # Run eigenvalue calculation and store results
            self.solver.computeEigenvalue(max_iters=self.max_iters)
            self.keffs.append(self.solver.getKeff())


    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Write a string with the results"""
        outstr = ''

        # Write out the mgxs file name eigenvalues from each simulation
        for file, keff in zip(self.mgxs_files, self.keffs):
            outstr += 'File: {0}\tkeff: {1:12.5E}\n'.format(file, keff)

        return outstr


if __name__ == '__main__':
    harness = LoadMGXSTestHarness()
    harness.main()
