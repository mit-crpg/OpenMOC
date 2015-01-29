## test homogeneous-one-group

## reliably works + passes

import numpy
from openmoc import *
import openmoc.log as log
from openmoc.options import Options
import unittest
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_directory[:-32])

from regression_test_runner import *
from benchmark_value_dict import *

def setup(sysargs):

    ###############################################################################
    #######################   Main Simulation Parameters   ########################
    ###############################################################################

    options = Options()

    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()

    log.set_log_level('NORMAL')

    log.py_printf('TITLE', 'Simulating a two group homogeneous infinite medium...')
    log.py_printf('HEADER', 'The reference keff = 1.72...')


    ###############################################################################
    ###########################   Creating Materials   ############################
    ###############################################################################

    log.py_printf('NORMAL', 'Creating materials...')

    infinite_medium = Material(name='2-group infinite medium')
    infinite_medium.setNumEnergyGroups(2)
    infinite_medium.setSigmaA(numpy.array([0.0038, 0.184]))
    infinite_medium.setSigmaF(numpy.array([0.000625, 0.135416667]))
    infinite_medium.setNuSigmaF(numpy.array([0.0015, 0.325]))
    infinite_medium.setSigmaS(numpy.array([0.1, 0.117, 0.0, 1.42]))
    infinite_medium.setChi(numpy.array([1.0, 0.0]))
    infinite_medium.setSigmaT(numpy.array([0.2208, 1.604]))


    ###############################################################################
    ###########################   Creating Surfaces   #############################
    ###############################################################################

    log.py_printf('NORMAL', 'Creating surfaces...')

    left = XPlane(x=-100.0, name='left')
    right = XPlane(x=100.0, name='right')
    top = YPlane(y=100.0, name='top')
    bottom = YPlane(y=-100.0, name='bottom')

    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)


    ###############################################################################
    #############################   Creating Cells   ##############################
    ###############################################################################

    log.py_printf('NORMAL', 'Creating cells...')

    cell = CellBasic()
    cell.setMaterial(infinite_medium)
    cell.addSurface(halfspace=+1, surface=left)
    cell.addSurface(halfspace=-1, surface=right)
    cell.addSurface(halfspace=+1, surface=bottom)
    cell.addSurface(halfspace=-1, surface=top)


    ###############################################################################
    #                            Creating Universes
    ###############################################################################

    log.py_printf('NORMAL', 'Creating universes...')

    root_universe = Universe(name='root universe')
    root_universe.addCell(cell)


    ###############################################################################
    ##########################   Creating the Geometry   ##########################
    ###############################################################################

    log.py_printf('NORMAL', 'Creating geometry...')

    geometry = Geometry()
    geometry.setRootUniverse(root_universe)
    geometry.initializeFlatSourceRegions()


    ###############################################################################
    ########################   Creating the TrackGenerator   ######################
    ###############################################################################

    log.py_printf('NORMAL', 'Initializing the track generator...')

    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.setNumThreads(num_threads)
    track_generator.generateTracks()


    ###############################################################################
    ###########################   Running a Simulation   ##########################
    ###############################################################################

    solver = CPUSolver(geometry, track_generator)
    solver.setNumThreads(num_threads)
    solver.setSourceConvergenceThreshold(tolerance)
    solver.convergeSource(max_iters)

    # return Keff
    return solver.getKeff()

# assign values for use in test case instance
test_type = 'Keff'
benchmark = 'homogeneous_two_groups'
benchmark_value = benchmark_val
print benchmark_value, 'is the Keff found for', benchmark
error_margin = 0.0001
filename = 'homogeneous-two-groups.py'
#setup_func = setup_homogeneous_two_groups

test_H2G = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads='DEFAULT')
test_H2G_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads=1)

test_list = [(test_H2G, __name__), (test_H2G_1t, __name__)]

if __name__ == '__main__':

    output = open("H1G_failed_results.txt", "a") # create output file in case of failures
    H2G_test_suite = regression_test_suite(test_list, output)
    H2G_test_suite.run_tests()
