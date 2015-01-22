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

output = open("H1G_failed_results.txt", "a") # create output file in case of failures

def setup_homogeneous_two_groups(sysargs):

    ## This setup function is named by adding 'setup_' to the beginning of the benchmark name
    ## for this test -- see bottom of file. self.benchmark for this test is
    ## 'homogeneous_two_group', so the function is named as seen.

    ## All setup functions have unique names so that when they are called by the
    ## regression test runner, only one function has that name. 
    
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
    infinite_medium.setNumEnergyGroups(2)
    infinite_medium.setSigmaA(numpy.array([0.0038, 0.184]))
    infinite_medium.setSigmaF(numpy.array([0.000625, 0.135416667]))
    infinite_medium.setNuSigmaF(numpy.array([0.0015, 0.325]))
    infinite_medium.setSigmaS(numpy.array([0.1, 0.117, 0.0, 1.42]))
    infinite_medium.setChi(numpy.array([1.0, 0.0]))
    infinite_medium.setSigmaT(numpy.array([0.2208, 1.604]))

    # surfaces
    circle = Circle(x=0.0, y=0.0, radius=50.0)
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
    cells.append(CellBasic(universe=1, material=1))
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
benchmark = 'homogeneous_two_groups'
benchmark_value = 1.323158836364746 # 1.723158836364746
error_margin = 0.005
filename = 'homogeneous-two-groups.py'
setup_func = setup_homogeneous_two_groups

test_H2G = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads='DEFAULT')
test_H2G_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, setup_func, num_threads=1)

test_list = [test_H2G, test_H2G_1t]

if __name__ == '__main__':
    # load case into suite (handles output file more completely)
    H2G_test_suite = regression_test_suite(test_list, output)
    H2G_test_suite.run_tests()
