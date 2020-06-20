#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
import openmoc

try:
    import openmc.openmoc_compatible
    import openmc.mgxs
except:
    print("OpenMC could not be imported, it's required for loading MGXS files")


class ComputeMaterialSPHTestHarness(TestHarness):
    """Computing SPH factors from a material-based OpenMC MGXS library using
       eigenvalue calculations. """

    def __init__(self):
        super(ComputeMaterialSPHTestHarness, self).__init__()
        self.azim_spacing = 0.12

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _create_solver(self):
        pass

    def _run_openmoc(self):
        """Run an SPH calculation."""

        # Initialize 16-group OpenMC multi-group XS library for a pin cell
        mgxs_lib = openmc.mgxs.Library.load_from_file(filename='mgxs',
             directory='../test_superhomogenization_material')

        # Compute SPH factors
        sph, sph_mgxs_lib, sph_indices = \
            openmoc.materialize.compute_sph_factors(
                mgxs_lib, azim_spacing=self.azim_spacing,
                num_azim=self.num_azim, num_threads=self.num_threads,
                sph_mode="eigenvalue", normalization="flux")

        # Create an OpenMOC Geometry from the OpenMOC Geometry
        openmoc_geometry = \
            openmc.openmoc_compatible.get_openmoc_geometry(sph_mgxs_lib.geometry)

        # Load cross section data
        openmoc_materials = \
            openmoc.materialize.load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

        # Initialize FSRs
        openmoc_geometry.initializeFlatSourceRegions()

        # Initialize an OpenMOC TrackGenerator
        track_generator = openmoc.TrackGenerator(
            openmoc_geometry, self.num_azim, self.azim_spacing)
        track_generator.generateTracks()

        # Initialize a new OpenMOC Solver
        self.solver = openmoc.CPUSolver(track_generator)
        self.solver.setConvergenceThreshold(self.tolerance)
        self.solver.setNumThreads(self.num_threads)

        # Compute an eigenvalue with the SPH library
        self.solver.computeEigenvalue()

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver"""
        return super(ComputeMaterialSPHTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = ComputeMaterialSPHTestHarness()
    harness.main()
