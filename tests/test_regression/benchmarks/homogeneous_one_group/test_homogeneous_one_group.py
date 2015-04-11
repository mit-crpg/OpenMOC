import numpy
from openmoc import *
import openmoc.log as log
from openmoc.options import Options
import unittest
import sys
import os

# these lines make sure the test_regression directory is in the path so regression_test_runner can be imported
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

    log.set_log_level('ERROR')
    

    ###############################################################################
    #######################   Main Simulation Parameters   ########################
    ###############################################################################


    infinite_medium = Material(name='1-group infinite medium')
    infinite_medium.setNumEnergyGroups(1)
    infinite_medium.setSigmaA(numpy.array([0.069389522]))
    infinite_medium.setSigmaF(numpy.array([0.0414198575]))
    infinite_medium.setNuSigmaF(numpy.array([0.0994076580]))
    infinite_medium.setSigmaS(numpy.array([0.383259177]))
    infinite_medium.setChi(numpy.array([1.0]))
    infinite_medium.setSigmaT(numpy.array([0.452648699]))


    ###############################################################################
    ###########################   Creating Surfaces   #############################
    ###############################################################################

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

    cell = CellBasic()
    cell.setMaterial(infinite_medium)
    cell.addSurface(halfspace=+1, surface=left)
    cell.addSurface(halfspace=-1, surface=right)
    cell.addSurface(halfspace=+1, surface=bottom)
    cell.addSurface(halfspace=-1, surface=top)


    ###############################################################################
    #                            Creating Universes
    ###############################################################################

    root_universe = Universe(name='root universe')
    root_universe.addCell(cell)

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

# assign values for use in test case instance
test_type = 'Keff'
benchmark = 'homogeneous_one_group'
benchmark_value = benchmark_val
print benchmark_value, 'is the Keff found for', benchmark
error_margin = 0.0001
filename = 'homogeneous-one-group.py'
#setup_func = setup_homogeneous_one_group

test_H1G = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, num_threads = 'DEFAULT')
test_H1G_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, num_threads = 1)

test_list = [(test_H1G, __name__), (test_H1G_1t, __name__)]

if __name__ == '__main__':
    
    output = open("H1G_failed_results.txt", "a") # create output file in case of failures
    H1G_test_suite = regression_test_suite(test_list, output)
    H1G_test_suite.run_tests()
    
