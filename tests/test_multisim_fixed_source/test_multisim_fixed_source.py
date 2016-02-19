#!/usr/bin/env python

import os
import sys
import hashlib
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import MultiSimTestHarness
from input_set import HomInfMedInput
import openmoc
import openmoc.process


class MultiSimFixedSourceTestHarness(MultiSimTestHarness):
    """A multi-simulation fixed source calculation for a cube of water with
    7-group cross sections and a fixed box source."""

    def __init__(self):
        super(MultiSimFixedSourceTestHarness, self).__init__()
        self.input_set = HomInfMedInput()
        self.res_type = openmoc.TOTAL_SOURCE
        self.solution_type = 'source'
        self.source_cell = None
        self.fluxes = []

    def _create_geometry(self):
        """Put a box source in the cube."""

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
        super(MultiSimFixedSourceTestHarness, self)._create_solver()
        self.solver.setFixedSourceByCell(self.source_cell, 1, 1.0)

    def _run_openmoc(self):
        """Run multiple OpenMOC fixed source calculations."""

        for i in range(self.num_simulations):
            super(MultiSimTestHarness, self)._run_openmoc()            
            self.num_iters.append(self.solver.getNumIterations())
            self.fluxes.append(openmoc.process.get_scalar_fluxes(self.solver))

    def _get_results(self, num_iterations=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=True):
        """Return eigenvalues from each simulation into a string."""

        # Write out the iteration count and fluxes from each simulation
        outstr = ''
        for i, num_iters in enumerate(self.num_iters):
            outstr += 'Iters: {0}\n'.format(num_iters)

            # Create a list of the floating point flux values
            fluxes = \
                ['{0:12.6E}'.format(flux) for flux in self.fluxes[i].ravel()]

            # Add the fluxes to the output string
            outstr += 'fluxes:\n'
            outstr += '\n'.join(fluxes) + '\n'

        # Hash the results if necessary.
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


if __name__ == '__main__':
    harness = MultiSimFixedSourceTestHarness()
    harness.main()
