#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc


class OTFTransportTestHarness(TestHarness):
    """An eigenvalue calculation for a 17x17 lattice with 7-group C5G7
    cross section data."""

    def __init__(self):
        super(OTFTransportTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(num_dimensions=3)
        self.num_polar = 4
        self.azim_spacing = 0.4
        self.z_spacing = 1.2
        self.tolerance = 1E-3

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(OTFTransportTestHarness, self)._create_geometry()

        # Change boundary condition on one side to vacuum
        cell_id = list(self.input_set.geometry.getRootUniverse().getCells())[0]
        surface_ids = list(self.input_set.geometry.getRootUniverse().getCells()[
             cell_id].getSurfaces())[1:5]
        del surface_ids[2]
        for surface_id in surface_ids:
            self.input_set.geometry.getRootUniverse().getCells()[
                 cell_id].getSurfaces()[surface_id].getSurface().setBoundaryType(
                 openmoc.VACUUM)

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(1.0)
        cmfd.setSORRelaxationFactor(1.5)
        cmfd.setLatticeStructure(4,4,4)
        cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
        cmfd.setKNearest(3)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setSegmentFormation(openmoc.OTF_STACKS)

    def _create_solver(self):
        """Instantiate a CPUSolver."""
        super(OTFTransportTestHarness, self)._create_solver()
        self.solver.setOTFTransport()

    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Need to use more than 1 thread, and lose reproducibility of FSR
        # numbering, in order to have temporary tracks and segments array of
        # the correct size for the multi-threaded solver.
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()

    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return hash as a string."""
        return super(OTFTransportTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)


if __name__ == '__main__':
    harness = OTFTransportTestHarness()
    harness.main()
