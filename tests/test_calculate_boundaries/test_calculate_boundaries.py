#!/usr/bin/env python

import os
import sys
import glob
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
import openmoc
from openmoc.log import py_printf


class CalculateBoundariesTestHarness(TestHarness):
    """A unit test of OpenMOC's Universe::calculateBoundaries routine."""

    def __init__(self):
        super(CalculateBoundariesTestHarness, self).__init__()

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _create_solver(self):
        pass

    def _run_openmoc(self):

        openmoc.set_log_level('NORMAL')

        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.XPlane(-40.)
        s2 = openmoc.XPlane(-50.)
        s3 = openmoc.XPlane(-60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(+1,s1)
        c2.addSurface(+1,s2)
        c2.addSurface(+1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMinXBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MinX: %f', u.getMinX())
        py_printf('NORMAL', 'MinXBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')

        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.YPlane(-40.)
        s2 = openmoc.YPlane(-50.)
        s3 = openmoc.YPlane(-60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(+1,s1)
        c2.addSurface(+1,s2)
        c2.addSurface(+1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMinYBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MinY: %f', u.getMinY())
        py_printf('NORMAL', 'MinYBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')


        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.ZPlane(-40.)
        s2 = openmoc.ZPlane(-50.)
        s3 = openmoc.ZPlane(-60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(+1,s1)
        c2.addSurface(+1,s2)
        c2.addSurface(+1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMinZBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MinZ: %f', u.getMinZ())
        py_printf('NORMAL', 'MinZBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')


        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.XPlane(40.)
        s2 = openmoc.XPlane(50.)
        s3 = openmoc.XPlane(60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(-1,s1)
        c2.addSurface(-1,s2)
        c2.addSurface(-1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMaxXBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MaxX: %f', u.getMaxX())
        py_printf('NORMAL', 'MaxXBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')

        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.YPlane(40.)
        s2 = openmoc.YPlane(50.)
        s3 = openmoc.YPlane(60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(-1,s1)
        c2.addSurface(-1,s2)
        c2.addSurface(-1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMaxYBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MaxY: %f', u.getMaxY())
        py_printf('NORMAL', 'MaxYBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')   

        c1 = openmoc.Cell()
        c2 = openmoc.Cell()

        s1 = openmoc.ZPlane(40.)
        s2 = openmoc.ZPlane(50.)
        s3 = openmoc.ZPlane(60.)
        s1.setBoundaryType(openmoc.VACUUM)
        s2.setBoundaryType(openmoc.PERIODIC)
        s3.setBoundaryType(openmoc.REFLECTIVE)

        c1.addSurface(-1,s1)
        c2.addSurface(-1,s2)
        c2.addSurface(-1,s3)

        u = openmoc.Universe()
        u.addCell(c1)
        u.addCell(c2)
        boundary = u.getMaxZBoundaryType()
        if boundary == 0:
            boundary = 'VACUUM'
        elif boundary == 1:
            boundary = 'REFLECTIVE'
        elif boundary == 2:
            boundary = 'PERIODIC'
        else:
            boundary = 'BOUNDARY_NONE'

        py_printf('NORMAL', 'MaxZ: %f', u.getMaxZ())
        py_printf('NORMAL', 'MaxZBoundaryType: %s', boundary)
        py_printf('SEPARATOR','')

        sW = openmoc.XPlane(10)
        sE = openmoc.XPlane(20)
        sS = openmoc.YPlane(30)
        sN = openmoc.YPlane(40)
        sB = openmoc.ZPlane(50)
        sT = openmoc.ZPlane(60)
        sW.setBoundaryType(openmoc.VACUUM)
        sE.setBoundaryType(openmoc.REFLECTIVE)
        sS.setBoundaryType(openmoc.VACUUM)
        sN.setBoundaryType(openmoc.REFLECTIVE)
        sB.setBoundaryType(openmoc.PERIODIC)
        sT.setBoundaryType(openmoc.REFLECTIVE)

        sX_mid = openmoc.XPlane(15)
        sY_mid = openmoc.YPlane(35)
        sZ_mid = openmoc.ZPlane(55)
        sX_mid.setBoundaryType(openmoc.BOUNDARY_NONE)
        sY_mid.setBoundaryType(openmoc.BOUNDARY_NONE)
        sZ_mid.setBoundaryType(openmoc.BOUNDARY_NONE)

        cell = openmoc.Cell()
        cell.addSurface(+1,sW)
        cell.addSurface(-1,sE)
        cell.addSurface(+1,sS)
        cell.addSurface(-1,sN)
        cell.addSurface(+1,sB)
        cell.addSurface(-1,sT)

        cell.addSurface(+1,sX_mid)
        cell.addSurface(-1,sX_mid)
        cell.addSurface(+1,sY_mid)
        cell.addSurface(-1,sY_mid)
        cell.addSurface(+1,sZ_mid)
        cell.addSurface(-1,sZ_mid)

        univ = openmoc.Universe()
        univ.addCell(cell)

        py_printf('NORMAL', 'MinX: %f', univ.getMinX())
        py_printf('NORMAL', 'MinXBoundaryType: %s', univ.getMinXBoundaryType())
        py_printf('NORMAL', 'MinY: %f', univ.getMinY())
        py_printf('NORMAL', 'MinYBoundaryType: %s', univ.getMinYBoundaryType())
        py_printf('NORMAL', 'MinZ: %f', univ.getMinZ())
        py_printf('NORMAL', 'MinZBoundaryType: %s', univ.getMinZBoundaryType())
        py_printf('NORMAL', 'MaxX: %f', univ.getMaxX())
        py_printf('NORMAL', 'MaxXBoundaryType: %s', univ.getMaxXBoundaryType())
        py_printf('NORMAL', 'MaxY: %f', univ.getMaxY())
        py_printf('NORMAL', 'MaxYBoundaryType: %s', univ.getMaxYBoundaryType())
        py_printf('NORMAL', 'MaxZ: %f', univ.getMaxZ())
        py_printf('NORMAL', 'MaxZBoundaryType: %s', univ.getMaxZBoundaryType())

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
    harness = CalculateBoundariesTestHarness()
    harness.main()
