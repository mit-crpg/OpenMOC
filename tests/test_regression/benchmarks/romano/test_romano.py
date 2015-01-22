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

output = open("romano_failed_results.txt", "a") # create output file in case of failures


def setup_romano(sysargs):

    sys.argv = sysargs

    ## parameters
    
    options = Options()

    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()

    log.set_log_level('ERROR')

    ## create materials

    fuel = Material(1)
    moderator = Material(2)

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

    ## create surfaces
    
    circle = Circle(x=0.0, y=0.0, radius=0.4)
    left = XPlane(x=-0.635)
    right = XPlane(x=0.635)
    top = YPlane(y=0.635)
    bottom = YPlane(y=-0.635)

    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)

    ## create cells

    cells = []
    cells.append(CellBasic(universe=1, material=1))
    cells.append(CellBasic(universe=1, material=2))
    cells.append(CellFill(universe=0, universe_fill=2))

    cells[0].addSurface(halfspace=-1, surface=circle)
    cells[1].addSurface(halfspace=+1, surface=circle)
    cells[2].addSurface(halfspace=+1, surface=left)
    cells[2].addSurface(halfspace=-1, surface=right)
    cells[2].addSurface(halfspace=+1, surface=bottom)
    cells[2].addSurface(halfspace=-1, surface=top)

    ## create lattice

    lattice = Lattice(id=2, width_x=1.27, width_y=1.27)
    lattice.setLatticeCells([[1]])

    ## create geometry

    geometry = Geometry()
    geometry.addMaterial(fuel)
    geometry.addMaterial(moderator)
    geometry.addCell(cells[0])
    geometry.addCell(cells[1])
    geometry.addCell(cells[2])
    geometry.addLattice(lattice)

    geometry.initializeFlatSourceRegions()

    ## create TrackGenerator

    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.generateTracks()

    ## run simulation

    solver = CPUSolver(geometry, track_generator)
    solver.setNumThreads(num_threads)
    solver.setSourceConvergenceThreshold(tolerance)
    solver.convergeSource(max_iters)

    return solver.getKeff()


#### assign values for use in creating the test case object
##test_type = 'Keff'
##benchmark = 'romano'
##benchmark_value = 1.2826626300811768
##Keff = general_romano_setup(['romano.py'])
##Keff_1t = general_romano_setup(['romano.py', '--num-omp-threads', '1'])
##error_margin = 0.0005
##
##test_romano = regression_test_case(test_type, benchmark, benchmark_value, Keff, Keff_1t, error_margin)
##
##if __name__ == '__main__':
##    # load case into suite (handles output file more completely)
##    romano_test_suite = regression_test_suite([test_romano], output)
##    romano_test_suite.run_tests()

# assign values for use in test case instance
test_type = 'Keff'
benchmark = 'romano'
benchmark_value = 1.1826626300811768 # 1.2826626300811768
error_margin = 0.005
filename = 'romano.py'
setup_func = setup_romano

test_romano = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads='DEFAULT')
test_romano_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads=1)

test_list = [test_romano, test_romano_1t]

if __name__ == '__main__':
    
    romano_test_suite = regression_test_suite(test_list, output)
    romano_test_suite.run_tests()


