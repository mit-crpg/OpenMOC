#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc
import math


class FindCellTestHarness(TestHarness):
    """Finding cells in a complex geometry"""

    def __init__(self):
        super(FindCellTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self.dimensions = 2

    def _create_geometry(self):
        """Inspired from SimpleLattice"""

        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymin = openmoc.YPlane(y=-2.0, name='ymin')
        ymax = openmoc.YPlane(y=+2.0, name='ymax')
        boundaries = [xmin, xmax, ymin, ymax]
        if (self.dimensions == 3):
            zmin = openmoc.ZPlane(z=-5.0, name='zmin')
            zmax = openmoc.ZPlane(z=+5.0, name='zmax')
            boundaries = [xmin, xmax, ymin, ymax, zmin, zmax]

        large_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                            radius=0.4, name='large pin')
        medium_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                             radius=0.3, name='medium pin')
        small_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                            radius=0.2, name='small pin')
        xplane = openmoc.XPlane(x=0, name='mid-cell')
        yplane = openmoc.YPlane(y=0, name='mid-cell')

        for boundary in boundaries: boundary.setBoundaryType(openmoc.REFLECTIVE)
        if (self.dimensions == 3):
            boundaries[-1].setBoundaryType(openmoc.VACUUM)

        large_fuel = openmoc.Cell(name='large pin fuel')
        large_fuel.setNumRings(3)
        large_fuel.setNumSectors(8)
        large_fuel.setFill(self.materials['UO2'])
        large_fuel.addSurface(halfspace=-1, surface=large_zcylinder)

        large_moderator = openmoc.Cell(name='large pin moderator')
        large_moderator.setNumSectors(8)
        large_moderator.setFill(self.materials['Water'])
        large_moderator.addSurface(halfspace=+1, surface=large_zcylinder)

        medium_fuel = openmoc.Cell(name='medium pin fuel')
        medium_fuel.setNumRings(3)
        medium_fuel.setNumSectors(8)
        medium_fuel.setFill(self.materials['UO2'])
        medium_fuel.addSurface(halfspace=-1, surface=medium_zcylinder)

        medium_moderator = openmoc.Cell(name='medium pin moderator')
        medium_moderator.setNumSectors(8)
        medium_moderator.setFill(self.materials['Water'])
        medium_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)

        small_fuel = openmoc.Cell(name='small pin fuel')
        small_fuel.setNumRings(3) # test may fail if fsrs are initialized
        small_fuel.setNumSectors(8)
        small_fuel.setFill(self.materials['UO2'])
        small_fuel.addSurface(halfspace=-1, surface=small_zcylinder)

        small_clad_q1 = openmoc.Cell(name='small pin cladding quarter')
        small_clad_q1.setFill(self.materials['Clad'])
        small_clad_q1.addSurface(halfspace=+1, surface=small_zcylinder)
        small_clad_q1.addSurface(halfspace=-1, surface=medium_zcylinder)
        small_clad_q1.addSurface(halfspace=+1, surface=xplane)
        small_clad_q1.addSurface(halfspace=-1, surface=yplane)

        small_moderator = openmoc.Cell(name='small pin moderator')
        small_moderator.setNumSectors(8)
        small_moderator.setFill(self.materials['Water'])
        small_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)

        lattice_cell = openmoc.Cell(name='lattice cell')

        root_cell = openmoc.Cell(name='root cell')
        root_cell.addSurface(halfspace=+1, surface=boundaries[0])
        root_cell.addSurface(halfspace=-1, surface=boundaries[1])
        root_cell.addSurface(halfspace=+1, surface=boundaries[2])
        root_cell.addSurface(halfspace=-1, surface=boundaries[3])
        if (self.dimensions == 3):
            root_cell.addSurface(halfspace=+1, surface=boundaries[4])
            root_cell.addSurface(halfspace=-1, surface=boundaries[5])

        pin1 = openmoc.Universe(name='large pin cell')
        pin2 = openmoc.Universe(name='medium pin cell')
        pin3 = openmoc.Universe(name='small pin cell')
        assembly = openmoc.Universe(name='2x2 lattice')
        root_universe = openmoc.Universe(name='root universe')

        pin1.addCell(large_fuel)
        pin1.addCell(large_moderator)
        pin2.addCell(medium_fuel)
        pin2.addCell(medium_moderator)
        pin3.addCell(small_fuel)
        pin3.addCell(small_clad_q1)
        pin3.addCell(small_moderator)
        assembly.addCell(lattice_cell)
        root_universe.addCell(root_cell)

        # 2x2 assembly
        lattice = openmoc.Lattice(name='2x2 lattice')
        lattice.setWidth(width_x=1.0, width_y=1.0)
        lattice.setUniverses([[[pin1, pin2], [pin1, pin3]]])
        lattice_cell.setFill(lattice)

        # Rotate a lattice to test rotation of a fill cell
        assembly_c = openmoc.Cell(name='rotated cell containing a lattice')
        assembly_c.setFill(assembly)
        assembly_c.setRotation((0,0,90), 'degrees')
        assembly_r = openmoc.Universe(name='rotated lattice')
        assembly_r.addCell(assembly_c)

        # Translate a lattice to test translation of a fill cell
        assembly_c = openmoc.Cell(name='rotated cell containing a lattice')
        assembly_c.setFill(assembly)
        assembly_c.setTranslation((0.05,0,0))
        assembly_t = openmoc.Universe(name='translated lattice')
        assembly_t.addCell(assembly_c)

        # Rotate a lattice to test rotation of a fill cell
        assembly_c = openmoc.Cell(name='translated cell containing a lattice')
        assembly_c.setFill(assembly)
        assembly_c.setRotation((0,0,-90), 'degrees')
        assembly_c.setTranslation((-0.05,0,0))
        assembly_rt = openmoc.Universe(name='translate and rotated lattice')
        assembly_rt.addCell(assembly_c)

        # 2x2 core
        core = openmoc.Lattice(name='2x2 core')
        core.setWidth(width_x=2.0, width_y=2.0)
        core.setUniverses([[[assembly, assembly_rt], [assembly_t, assembly_r]]])
        root_cell.setFill(core)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

    def _create_solver(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _run_openmoc(self):
        pass

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return as a string."""
        outstr = ''

        # Short hands
        geometry = self.geometry
        root = geometry.getRootUniverse()
        rad = 0.2
        rad_f = 0.1
        pi = math.pi

        # Find cells in the complex geometry
        # Regular cell in lattice (top left lattice)
        coords = openmoc.LocalCoords(-0.5+rad, 0.5-rad, 0)
        coords.setUniverse(root)
        cell = geometry.findCellContainingCoords(coords)
        outstr += cell.toString() + '\n'

        # Translated cell (bottom left lattice)
        coords = openmoc.LocalCoords(-0.5+rad+0.05, -1.5-rad, 0)
        coords.setUniverse(root)
        cell = geometry.findCellContainingCoords(coords)
        outstr += cell.toString() + '\n'

        # Rotated cell (bottom right lattice)
        coords = openmoc.LocalCoords(1.5+rad, -0.5+rad, 0)
        coords.setUniverse(root)
        cell = geometry.findCellContainingCoords(coords)
        outstr += cell.toString() + '\n'

        # Find the next cell in a direction in the complex geometry. First, the
        # LocalCoords have to be populated with findCellContainingCoords
        # Regular cell in lattice
        coords = openmoc.LocalCoords(-0.5+rad_f, 0.5-rad_f, 0)
        coords.setUniverse(root)
        geometry.findCellContainingCoords(coords)
        cell = geometry.findNextCell(coords, -pi/4, pi/4)
        outstr += cell.toString() + '\n'

        # Translated cell
        coords = openmoc.LocalCoords(-0.5+rad_f+0.05, -1.5-rad_f, 0)
        coords.setUniverse(root)
        geometry.findCellContainingCoords(coords)
        cell = geometry.findNextCell(coords, -pi/4, pi/2)
        outstr += cell.toString() + '\n'

        # Rotated cell
        coords = openmoc.LocalCoords(1.5+rad_f, -0.5+rad_f, 0)
        coords.setUniverse(root)
        geometry.findCellContainingCoords(coords)
        cell = geometry.findNextCell(coords, +pi/4, 3*pi/4)
        outstr += cell.toString() + '\n'

        return outstr


if __name__ == '__main__':
    harness = FindCellTestHarness()
    harness.main()
