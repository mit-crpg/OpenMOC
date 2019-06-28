#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import PinCellInput
import openmoc


class PinCellTestHarness(TestHarness):
    """An eigenmodes calculation in a pin cell with 7-group C5G7 data."""
    def __init__(self):
        super(PinCellTestHarness, self).__init__()
        self.input_set = PinCellInput()

    def _create_geometry(self):
        """Instantiate materials and a pin cell Geometry."""

        self.input_set.create_materials()

        zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymin = openmoc.YPlane(y=-2.0, name='ymin')
        ymax = openmoc.YPlane(y=+2.0, name='ymax')

        xmin.setBoundaryType(openmoc.VACUUM)
        xmax.setBoundaryType(openmoc.VACUUM)
        ymin.setBoundaryType(openmoc.VACUUM)
        ymax.setBoundaryType(openmoc.VACUUM)

        fuel = openmoc.Cell(name='fuel')
        fuel.setFill(self.input_set.materials['UO2'])
        fuel.addSurface(halfspace=-1, surface=zcylinder)

        moderator = openmoc.Cell(name='moderator')
        moderator.setFill(self.input_set.materials['Water'])
        moderator.addSurface(halfspace=+1, surface=zcylinder)
        moderator.addSurface(halfspace=+1, surface=xmin)
        moderator.addSurface(halfspace=-1, surface=xmax)
        moderator.addSurface(halfspace=+1, surface=ymin)
        moderator.addSurface(halfspace=-1, surface=ymax)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(fuel)
        root_universe.addCell(moderator)

        self.input_set.geometry = openmoc.Geometry()
        self.input_set.geometry.setRootUniverse(root_universe)

    def _create_solver(self):
        """Instantiate a IRAMSolver."""
        self.solver = openmoc.CPUSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)

        # Initialize IRAMSolver to perform forward eigenmode calculation
        self.solver = openmoc.krylov.IRAMSolver(self.solver)

    def _run_openmoc(self):
        """Run an OpenMOC forward kyrlov calculation."""
        self.solver.computeEigenmodes(num_modes=1, solver_mode=openmoc.ADJOINT)

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return as a string."""

        # Round to 10th decimal to avoid floating point issues
        outstr = str([round(k, 10) for k in self.solver._eigenvalues])
        
        return outstr

if __name__ == '__main__':
    harness = PinCellTestHarness()
    harness.main()
