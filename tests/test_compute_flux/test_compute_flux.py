#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import HomInfMedInput
import openmoc


class ComputeFluxTestHarness(TestHarness):
    """A fixed source flux calculation in a cube of water with 7-group
    cross sections. This tests the Solver::computeFlux(...) method."""

    def __init__(self):
        super(ComputeFluxTestHarness, self).__init__()
        self.input_set = HomInfMedInput()
        self.res_type = openmoc.SCALAR_FLUX
        self.solution_type = 'flux'

    def _create_geometry(self):
        """Put a box source """

        self.input_set.create_materials()
        self.input_set.create_geometry()

        # Get the root Cell
        cells = self.input_set.geometry.getAllCells()
        for cell_id in cells:
            cell = cells[cell_id]
            if cell.getName() == 'root cell':
                root_cell = cell

        # Apply VACUUM BCs on all bounding surfaces
        surfaces = root_cell.getSurfaces()
        for surface_id in surfaces:
            surface = surfaces[surface_id]._surface
            surface.setBoundaryType(openmoc.VACUUM)

        # Replace fissionable infinite medium material with C5G7 water
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

        lattice = openmoc.castUniverseToLattice(root_cell.getFillUniverse())
        num_x = lattice.getNumX()
        num_y = lattice.getNumY()
        width_x = lattice.getWidthX()
        width_y = lattice.getWidthY()

        # Create cells filled with water to put in Lattice
        water_cell = openmoc.Cell(name='water')
        water_cell.setFill(self.materials['Water'])
        water_univ = openmoc.Universe(name='water')
        water_univ.addCell(water_cell)

        self.source_cell = openmoc.Cell(name='source')
        self.source_cell.setFill(self.materials['Water'])
        source_univ = openmoc.Universe(name='source')
        source_univ.addCell(self.source_cell)

        # Create 2D array of Universes in each lattice cell
        universes = [[[water_univ]*num_x for _ in range(num_y)]]

        # Place fixed source Universe at (x=0.5, y=0.5)
        source_x = 0.5
        source_y = 0.5
        lat_x = (root_cell.getMaxX() - source_x) / width_x
        lat_y = (root_cell.getMaxY() - source_y) / width_y
        universes[0][int(lat_x)][int(lat_y)] = source_univ

        # Create a new Lattice for the Universes
        lattice = openmoc.Lattice(name='{0}x{1} lattice'.format(num_x, num_y))
        lattice.setWidth(width_x=width_x, width_y=width_y)
        lattice.setUniverses(universes)
        root_cell.setFill(lattice)

    def _create_solver(self):
        """Instantiate a CPUSolver."""
        super(ComputeFluxTestHarness, self)._create_solver()
        self.solver.setFixedSourceByCell(self.source_cell, 1, 1.0)
        self.solver.setFixedSourceByCell(self.source_cell, 2, 0.5)
        self.solver.setFixedSourceByCell(self.source_cell, 3, 0.25)

        geometry = self.input_set.geometry
        track_generator = self.track_generator
        solver=self.solver


if __name__ == '__main__':
    harness = ComputeFluxTestHarness()
    harness.main()
