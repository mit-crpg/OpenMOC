#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import InputSet
import openmoc



class SimplerLatticeInput(InputSet):
    """A 4x4 pin cell lattice problem from sample-input/simple-lattice."""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_lattice(self):
    
        large_zcylinder = openmoc.ZCylinder(
            x=0.0, y=0.0, radius=0.4, name='large pin')
        medium_zcylinder = openmoc.ZCylinder(
            x=0.0, y=0.0, radius=0.3, name='medium pin')
        small_zcylinder = openmoc.ZCylinder(
            x=0.0, y=0.0, radius=0.2, name='small pin')
        x0_plane = openmoc.XPlane(x=0)
        y0_plane = openmoc.YPlane(y=0)
    
        large_fuel = openmoc.Cell(name='large pin fuel')
        large_fuel.setFill(self.materials['UO2'])
        large_fuel.addSurface(halfspace=-1, surface=large_zcylinder)
    
        large_moderator = openmoc.Cell(name='large pin moderator')
        large_moderator.setFill(self.materials['Water'])
        large_moderator.addSurface(halfspace=+1, surface=large_zcylinder)
        
        fuel_quads = [None]*4
        mod_quads = [None]*4
        for i in range(4):
            qname = " (quadrant {})".format(i + 1)
            fcell = openmoc.Cell(name='large fuel' + qname)
            fcell.setFill(self.materials['UO2'])
            fcell.addSurface(halfspace=-1, surface=large_zcylinder)
            mcell = openmoc.Cell(name='large mod' + qname)
            mcell.setFill(self.materials['Water'])
            mcell.addSurface(halfspace=+1, surface=large_zcylinder)
            # Add the appropriate quadrants
            xspace = 1 if i in (0, 3) else -1
            yspace = 1 if i in (0, 1) else -1
            fcell.addSurface(halfspace=xspace, surface=x0_plane)
            fcell.addSurface(halfspace=yspace, surface=y0_plane)
            mcell.addSurface(halfspace=xspace, surface=x0_plane)
            mcell.addSurface(halfspace=yspace, surface=y0_plane)
            fuel_quads[i] = fcell
            mod_quads[i] = mcell
    
        medium_fuel = openmoc.Cell(name='medium pin fuel')
        medium_fuel.setFill(self.materials['UO2'])
        medium_fuel.addSurface(halfspace=-1, surface=medium_zcylinder)
    
        medium_moderator = openmoc.Cell(name='medium pin moderator')
        medium_moderator.setFill(self.materials['Water'])
        medium_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)
    
        small_fuel = openmoc.Cell(name='small pin fuel')
        small_fuel.setFill(self.materials['UO2'])
        small_fuel.addSurface(halfspace=-1, surface=small_zcylinder)
    
        small_moderator = openmoc.Cell(name='small pin moderator')
        small_moderator.setFill(self.materials['Water'])
        small_moderator.addSurface(halfspace=+1, surface=small_zcylinder)
    
        pin1 = openmoc.Universe(name='large pin cell')
        pin2 = openmoc.Universe(name='medium pin cell')
        pin3 = openmoc.Universe(name='small pin cell')
        pin4 = openmoc.Universe(name='large pin cell (quadrants)')
       
        pin1.addCell(large_fuel)
        pin1.addCell(large_moderator)
        pin2.addCell(medium_fuel)
        pin2.addCell(medium_moderator)
        pin3.addCell(small_fuel)
        pin3.addCell(small_moderator)
        for c in fuel_quads + mod_quads:
            pin4.addCell(c)
        
        ###################
        # SUBDIVIDE PINS  #
        ###################
        lleft = [-0.5]*2
        uright = [0.5]*2
        DIV = (2, 2)
        div1 = openmoc.Subdivider(DIV, lleft, uright)
        pin1_div = div1.get_subdivided_universe(pin1)
        div4 = openmoc.Subdivider(DIV, lleft, uright)
        pin4_div = div4.get_subdivided_universe(pin4)
    
        # 2x2 assembly
        lattice = openmoc.Lattice(name='2x2 lattice')
        lattice.setWidth(width_x=1.0, width_y=1.0)
        lattice.setUniverses([[[pin4_div,     pin3],
                               [pin1,     pin1_div]]])
        
        return lattice
    
    def create_geometry(self):
        """Instantiate a 4x4 pin cell lattice Geometry."""
        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymin = openmoc.YPlane(y=-2.0, name='ymin')
        ymax = openmoc.YPlane(y=+2.0, name='ymax')
        boundaries = [xmin, xmax, ymin, ymax]
        if (self.dimensions == 3):
            zmin = openmoc.ZPlane(z=-2.0, name='zmin')
            zmax = openmoc.ZPlane(z=+2.0, name='zmax')
            boundaries = [xmin, xmax, ymin, ymax, zmin, zmax]
        for boundary in boundaries:
            boundary.setBoundaryType(openmoc.REFLECTIVE)
        if (self.dimensions == 3):
            boundaries[-1].setBoundaryType(openmoc.VACUUM)

        lat = self.create_lattice()

        root_cell = openmoc.Cell(name='root cell')
        root_cell.addSurface(halfspace=+1, surface=boundaries[0])
        root_cell.addSurface(halfspace=-1, surface=boundaries[1])
        root_cell.addSurface(halfspace=+1, surface=boundaries[2])
        root_cell.addSurface(halfspace=-1, surface=boundaries[3])
        if (self.dimensions == 3):
            root_cell.addSurface(halfspace=+1, surface=boundaries[4])
            root_cell.addSurface(halfspace=-1, surface=boundaries[5])

        lattice_cell = openmoc.Cell(name='lattice cell')
        assembly = openmoc.Universe(name='2x2 lattice')
        root_universe = openmoc.Universe(name='root universe')
        assembly.addCell(lattice_cell)
        root_universe.addCell(root_cell)
        lattice_cell.setFill(lat)
        # 2x2 core
        core = openmoc.Lattice(name='2x2 core')
        core.setWidth(width_x=2.0, width_y=2.0)
        core.setUniverses([[[assembly, assembly], [assembly, assembly]]])
        root_cell.setFill(core)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(SimplerLatticeInput, self).create_geometry()

class SubdividingLatticeTestHarness(TestHarness):
    """An eigenvalue calculation for a 4x4 lattice with 7-group C5G7
    cross section data."""

    def __init__(self):
        super(SubdividingLatticeTestHarness, self).__init__()
        self.input_set = SimplerLatticeInput(num_dimensions=2)

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return hash as a string."""
        return super(SubdividingLatticeTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)


if __name__ == '__main__':
    harness = SubdividingLatticeTestHarness()
    harness.main()
