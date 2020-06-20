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


class ComputeCellSPHTestHarness(TestHarness):
    """Computing SPH factors from a cell-based OpenMC MGXS library. Factors
       are used in the fuel and water in this test case."""

    def __init__(self):
        super(ComputeCellSPHTestHarness, self).__init__()
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

        # Initialize 16-group OpenMC multi-group cross section library for a pin cell
        mgxs_lib = openmc.mgxs.Library.load_from_file(filename='mgxs', directory='.')

        # Create an OpenMOC Geometry from the OpenMOC Geometry
        openmoc_geometry = \
            openmc.openmoc_compatible.get_openmoc_geometry(mgxs_lib.geometry)

        # Load cross section data
        openmoc_materials = \
            openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

        # Create list of cells (fuel+water) which should have sph factors
        sph_domains = []
        for cell in openmoc_geometry.getAllMaterialCells().values():
            if "fuel" in cell.getName() or "water" in cell.getName():
                sph_domains.append(cell.getId())

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

        # Compute SPH factors
        sph, sph_indices = \
            openmoc.materialize.compute_sph_factors(
                mgxs_lib, azim_spacing=self.azim_spacing,
                num_azim=self.num_azim, num_threads=self.num_threads,
                solver=self.solver, track_generator=track_generator,
                geometry=openmoc_geometry, sph_domains=sph_domains,
                return_library=False)

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver"""
        return super(ComputeCellSPHTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = ComputeCellSPHTestHarness()
    harness.main()
