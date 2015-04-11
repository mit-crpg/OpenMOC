## Test Romano benchmark

import numpy
from openmoc import *
import openmoc.log as log
from openmoc.options import Options
import unittest
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_directory[:-17])

from regression_test_runner import *
from benchmark_value_dict import *

def setup_romano(sysargs):

    sys.argv = sysargs

    ###############################################################################
    #######################   Main Simulation Parameters   ########################
    ###############################################################################

    options = Options()

    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()

    log.set_log_level('ERROR')


    ###############################################################################
    ###########################   Creating Materials   ############################
    ###############################################################################

    fuel = Material(name='fuel')
    moderator = Material(name='moderator')

    fuel.setNumEnergyGroups(1)
    moderator.setNumEnergyGroups(1)

    fuel.setSigmaA(numpy.array([0.069389522]))
    fuel.setSigmaT(numpy.array([0.452648699]))
    fuel.setSigmaF(numpy.array([0.0414198575]))
    fuel.setNuSigmaF(numpy.array([0.0994076580]))
    fuel.setSigmaS(numpy.array([0.38259177]))
    fuel.setChi(numpy.array([1.0]))

    moderator.setSigmaA(numpy.array([0.003751099]))
    moderator.setSigmaT(numpy.array([0.841545641]))
    moderator.setSigmaF(numpy.array([0.0]))
    moderator.setNuSigmaF(numpy.array([0.0]))
    moderator.setSigmaS(numpy.array([0.837794542]))
    moderator.setChi(numpy.array([1.0]))


    ###############################################################################
    ###########################   Creating Surfaces   #############################
    ###############################################################################

    circle = Circle(x=0.0, y=0.0, radius=0.4)
    left = XPlane(x=-0.635)
    right = XPlane(x=0.635)
    top = YPlane(y=0.635)
    bottom = YPlane(y=-0.635)

    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)


    ###############################################################################
    #############################   Creating Cells   ##############################
    ###############################################################################

    fuel_cell = CellBasic(name='fuel')
    fuel_cell.setMaterial(fuel)
    fuel_cell.addSurface(halfspace=-1, surface=circle)

    moderator_cell = CellBasic(name='moderator')
    moderator_cell.setMaterial(moderator)
    moderator_cell.addSurface(halfspace=+1, surface=circle)
    moderator_cell.addSurface(halfspace=+1, surface=left)
    moderator_cell.addSurface(halfspace=-1, surface=right)
    moderator_cell.addSurface(halfspace=+1, surface=bottom)
    moderator_cell.addSurface(halfspace=-1, surface=top)


    ###############################################################################
    ###########################   Creating Universes   ############################
    ###############################################################################

    root_universe = Universe(name='root universe')
    root_universe.addCell(fuel_cell)
    root_universe.addCell(moderator_cell)


    ###############################################################################
    ##########################   Creating the Geometry   ##########################
    ###############################################################################

    geometry = Geometry()
    geometry.setRootUniverse(root_universe)
    geometry.initializeFlatSourceRegions()


    ###############################################################################
    ########################   Creating the TrackGenerator   ######################
    ###############################################################################

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

    return solver.getKeff()


test_type = 'Keff'
benchmark = 'romano'
benchmark_value = benchmark_val
print benchmark_value, 'is the Keff found for', benchmark
error_margin = 0.0001
filename = 'romano.py'
#setup_func = setup_romano

test_romano_case = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads='DEFAULT')
test_romano_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads=1)

test_list = [(test_romano_case, test_romano), (test_romano_1t, test_romano)]

if __name__ == '__main__':

    output = open("romano_failed_results.txt", "a") # create output file in case of failures
    romano_test_suite = regression_test_suite(test_list, output)
    romano_test_suite.run_tests()


