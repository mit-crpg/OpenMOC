#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from openmoc.log import py_printf
import os

class GeometryDumpTestHarness(TestHarness):
    """Test dumping a geometry to file."""

    def __init__(self):
        super(GeometryDumpTestHarness, self).__init__()


    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')


    def _create_geometry(self):
        self.create_materials()


    def _create_solver(self):
        pass


    def _create_trackgenerator(self):
        pass


    def _generate_tracks(self):
        pass


    def _run_openmoc(self):
        """Instantiate a complex region Geometry."""
        root_universe = openmoc.Universe(name='root universe')
        root_cell = openmoc.Cell(name='root cell')
        root_universe.addCell(root_cell)
        u1 = openmoc.Universe(name='universe 1')
        u2 = openmoc.Universe(name='universe 2 in c2')
        root_cell.setFill(u1)

        # Z-bounds at root level
        p1 = openmoc.ZPlane(z=-2.0, name='zmin')
        p1.setBoundaryType(openmoc.REFLECTIVE)
        p2 = openmoc.ZPlane(z=+2.0, name='zmax')
        p2.setBoundaryType(openmoc.INTERFACE)
        root_cell.addSurface(halfspace=+1, surface=p1)
        root_cell.addSurface(halfspace=-1, surface=p2)

        # Cells in the root cell
        c1 = openmoc.Cell(name='intersection cell')
        c2 = openmoc.Cell(name='union cell')
        c2a = openmoc.Cell(name='union cell 2')
        c3 = openmoc.Cell(name='unbound cell')
        u1.addCell(c1)
        u1.addCell(c2)
        u1.addCell(c3)
        # Setting the parent cell helps to find boundaries (in Z here)
        c3.setParent(root_cell)

        # Cell c2a in c2, to further test arborescence
        c2.setFill(u2)
        u2.addCell(c2a)
        c2.setParent(root_cell)
        c2a.setParent(c2)  # to test transitivity

        # Bounds for cell 1 : intersection region cell
        p11 = openmoc.XPlane(x=-2.0, name='xmin')
        p11.setBoundaryType(openmoc.REFLECTIVE)
        p12 = openmoc.YPlane(y=-2.0, name='ymin')
        p12.setBoundaryType(openmoc.VACUUM)
        p13 = openmoc.ZCylinder(x=0,y=0,radius=2.5, name='cylinder')
        p13.setBoundaryType(openmoc.INTERFACE)

        # addSurface assumes an intersection, which is what we want here
        c1.addSurface(halfspace=+1, surface=p11)
        c1.addSurface(halfspace=+1, surface=p12)
        c1.addSurface(halfspace=-1, surface=p13)

        # Bounds for cell 2 : union region cell
        p22 = openmoc.ZCylinder(x=4,y=4,radius=2.5, name='cylinder')
        p22.setBoundaryType(openmoc.INTERFACE)
        p23 = openmoc.ZCylinder(x=-2,y=-8,radius=3, name='cylinder')
        p23.setBoundaryType(openmoc.VACUUM)

        # To have a union, we need to use addLogicalNode and addSurfaceInRegion
        c2.addLogicalNode(1)
        c2.addSurfaceInRegion(halfspace=-1, surface=p22)
        c2.addSurfaceInRegion(halfspace=-1, surface=p23)

        # Plane limits area even more
        c2a.addLogicalNode(1)
        c2a.addSurfaceInRegion(halfspace=+1, surface=p11)

        # Add a material to a cell to know the number of groups
        c1.setFill(self.materials['UO2'])

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        # Dump geometry
        self.geometry.dumpToFile("geometry_file.geo")

        # Get rid of the geometry
        self.geometry = None

        # Reload the geometry
        self.geometry = openmoc.Geometry()
        self.geometry.loadFromFile("geometry_file.geo")

        # Dump geometry again
        self.geometry.dumpToFile("geometry_file_second.geo")

        # NOTE : dumping the geometry before and after FSRs are initialized
        # leads to different results, as the rings and sectors create new cells.


    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Compare the two geometry files."""

        # Assert that both files are the same
        if (os.system("cmp geometry_file.geo geometry_file_second.geo") != 0):
            py_printf('ERROR', "Geometry files are not dumped and loaded "
                      "properly")

        # Python2 and 3 binaries differ, cannot compare the file to a reference.
        outstr = 'dummy'
        return outstr

if __name__ == '__main__':
    harness = GeometryDumpTestHarness()
    harness.main()
