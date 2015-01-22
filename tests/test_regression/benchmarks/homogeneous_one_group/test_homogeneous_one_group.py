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

def setup_homogeneous_one_group(sysargs):
    
    # run simulation w/ given command line arguments
    sys.argv = sysargs
    options = Options()
    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()
    log.set_log_level('ERROR')
    
    # materials   
    infinite_medium = Material(1)
    infinite_medium.setNumEnergyGroups(1)
    infinite_medium.setSigmaA(numpy.array([0.069389522]))
    infinite_medium.setSigmaF(numpy.array([0.0414198575]))
    infinite_medium.setNuSigmaF(numpy.array([0.0994076580]))
    infinite_medium.setSigmaS(numpy.array([0.383259177]))
    infinite_medium.setChi(numpy.array([1.0]))
    infinite_medium.setSigmaT(numpy.array([0.452648699]))

    # surfaces
    circle = Circle(x=0.0, y=0.0, radius=10.0)
    left = XPlane(x=-100.0)
    right = XPlane(x=100.0)
    top = YPlane(y=100.0)
    bottom = YPlane(y=-100.0)

    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)

    # cells
    cells = []
    cells.append(CellBasic(universe=1, material=1, rings=2, sectors=4))
    cells.append(CellBasic(universe=1, material=1))
    cells.append(CellFill(universe=0, universe_fill=2))

    cells[0].addSurface(halfspace=-1, surface=circle)
    cells[1].addSurface(halfspace=+1, surface=circle)
    cells[2].addSurface(halfspace=+1, surface=left)
    cells[2].addSurface(halfspace=-1, surface=right)
    cells[2].addSurface(halfspace=+1, surface=bottom)
    cells[2].addSurface(halfspace=-1, surface=top)

    # lattices
    lattice = Lattice(id=2, width_x=200.0, width_y=200.0)
    lattice.setLatticeCells([[1]])

    # geometry
    geometry = Geometry()
    geometry.addMaterial(infinite_medium)
    geometry.addCell(cells[0])
    geometry.addCell(cells[1])
    geometry.addCell(cells[2])
    geometry.addLattice(lattice)

    geometry.initializeFlatSourceRegions()

    # TrackGenerator
    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.generateTracks()

    # run simulation
    solver = CPUSolver(geometry, track_generator)
    solver.setNumThreads(num_threads)
    solver.setSourceConvergenceThreshold(tolerance)
    solver.convergeSource(max_iters)

    # return Keff
    return solver.getKeff()

# assign values for use in test case instance
test_type = 'Keff'
benchmark = 'homogeneous_one_group'
benchmark_value = benchmark_value_dictionary[(benchmark,test_type)]
error_margin = 0.0001
filename = 'homogeneous-one-group.py'
setup_func = setup_homogeneous_one_group

test_H1G = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads = 'DEFAULT')
test_H1G_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads = 1)

test_list = [test_H1G, test_H1G_1t]

if __name__ == '__main__':
    
    output = open("H1G_failed_results.txt", "a") # create output file in case of failures
    H1G_test_suite = regression_test_suite(test_list, output)
    H1G_test_suite.run_tests()
