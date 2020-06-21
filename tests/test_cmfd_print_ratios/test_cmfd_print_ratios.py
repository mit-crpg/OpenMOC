#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PwrAssemblyInput
import openmoc
from openmoc.log import py_printf


class CmfdPwrAssemblyTestHarness(TestHarness):
    """An eigenvalue calculation for a 17x17 lattice with 7-group C5G7
    cross section data."""

    def __init__(self):
        super(CmfdPwrAssemblyTestHarness, self).__init__()
        self.input_set = PwrAssemblyInput()
        self.spacing = 0.12
        self.max_iters = 10

        #FIXME
        self.num_threads = 1

    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        super(CmfdPwrAssemblyTestHarness, self)._create_geometry()

        # Change boundary condition on two opposite sides to periodic
        cell_id = list(self.input_set.geometry.getRootUniverse().getCells())[0]
        surface_id = list(self.input_set.geometry.getRootUniverse().getCells()[
             cell_id].getSurfaces())[2]
        self.input_set.geometry.getRootUniverse().getCells()[
             cell_id].getSurfaces()[surface_id].getSurface().setBoundaryType(
             openmoc.PERIODIC)
        surface_id = list(self.input_set.geometry.getRootUniverse().getCells()[
             cell_id].getSurfaces())[3]
        self.input_set.geometry.getRootUniverse().getCells()[
             cell_id].getSurfaces()[surface_id].getSurface().setBoundaryType(
             openmoc.PERIODIC)

        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(0.4)
        cmfd.setSORRelaxationFactor(1.6)
        cmfd.setLatticeStructure(17,17)
        cmfd.setGroupStructure([[1],[2],[3],[4],[5],[6,7]])

        # Print CMFD update ratios after every iteration
        cmfd.printProlongation()

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)

    def _get_results(self, num_iters=False, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):

        """Compare an update ratio file"""
        # Assert that both files are the same
        if (os.system("cmp pf_group_4_iter_8.txt pf_reference.txt") != 0):
            py_printf('ERROR', "Prolongation ratios are different")

        """Digest info in the solver and return hash as a string."""
        return super(CmfdPwrAssemblyTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

    def _cleanup(self):
        """Delete the update ratio files."""
        super(CmfdPwrAssemblyTestHarness, self)._cleanup()
        os.system("rm pf_group_*_iter_*.txt")


if __name__ == '__main__':
    harness = CmfdPwrAssemblyTestHarness()
    harness.main()
