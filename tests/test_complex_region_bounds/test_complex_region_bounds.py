#!/usr/bin/env python

import os
import sys
import glob
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
import openmoc
from openmoc.log import py_printf


class ComplexRegionBoundsTestHarness(TestHarness):
    """Tests complex cells min and max coordinate/boundary type getters."""

    def __init__(self):
        super(ComplexRegionBoundsTestHarness, self).__init__()

    def _run_openmoc(self):
        """Instantiate a complex region Geometry."""
#                                ----------------
#                              /                 \
#             --------------- c2 (p22,p23) - u2 - c2a (p11)
#           /               /
#  u_r <- c_r (p1,p2) - u1 - c1 (p11,p12,p13)
#           \             \
#            ------------- c3
        root_universe = openmoc.Universe(name='root universe')
        root_cell = openmoc.Cell(name='root cell')
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

        openmoc.set_log_level('NORMAL')
        for cell in [root_cell, c1, c2, c2a, c3]:
            py_printf('NORMAL', 'Cell: %s', cell.getName())
            py_printf('NORMAL', 'MinX: %f', cell.getMinX())
            py_printf('NORMAL', 'MinXBoundaryType: %s', cell.getMinXBoundaryType())
            py_printf('NORMAL', 'MinY: %f', cell.getMinY())
            py_printf('NORMAL', 'MinYBoundaryType: %s', cell.getMinYBoundaryType())
            py_printf('NORMAL', 'MinZ: %f', cell.getMinZ())
            py_printf('NORMAL', 'MinZBoundaryType: %s', cell.getMinZBoundaryType())
            py_printf('NORMAL', 'MaxX: %f', cell.getMaxX())
            py_printf('NORMAL', 'MaxXBoundaryType: %s', cell.getMaxXBoundaryType())
            py_printf('NORMAL', 'MaxY: %f', cell.getMaxY())
            py_printf('NORMAL', 'MaxYBoundaryType: %s', cell.getMaxYBoundaryType())
            py_printf('NORMAL', 'MaxZ: %f', cell.getMaxZ())
            py_printf('NORMAL', 'MaxZBoundaryType: %s', cell.getMaxZBoundaryType())
            py_printf('NORMAL', '')

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _create_solver(self):
        pass

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the log file and return as a string."""

        # Find the log filename with the time and date
        logfilename = glob.glob('log/openmoc-*')

        # Read the file into a list of strings for each line
        with open(logfilename[0], 'r') as myfile:
            lines = myfile.readlines()

        # Concatenate all strings in the file into a single string
        # Exclude the first line which is the time and date
        outstr = ''.join(lines[1:])
        return outstr

if __name__ == '__main__':
    harness = ComplexRegionBoundsTestHarness()
    harness.main()
