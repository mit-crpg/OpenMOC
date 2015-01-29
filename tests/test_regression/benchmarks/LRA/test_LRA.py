import os

from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options
import unittest
import sys
import time

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_directory[:-15])

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
    ###########################   Creating Materials   ############################
    ###############################################################################

    materials = materialize.materialize(current_directory+'/LRA-materials.py')

    ###############################################################################
    ###########################   Creating Surfaces   #############################
    ###############################################################################

    left = XPlane(x=-82.5)
    right = XPlane(x=82.5)
    bottom = YPlane(y=-82.5)
    top = YPlane(y=82.5)
    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(VACUUM)
    bottom.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(VACUUM)

    ###############################################################################
    ######################   Creating Cells and Universes   #######################
    ###############################################################################

    # Region 1
    region1_cell = CellBasic(name='region 1')
    region1_cell.setMaterial(materials['region_1'])
    region1 = Universe(name='region 1')
    region1.addCell(region1_cell)

    # Region 2
    region2_cell = CellBasic(name='region 2')
    region2_cell.setMaterial(materials['region_2'])
    region2 = Universe(name='region 2')
    region2.addCell(region2_cell)

    # Region 3
    region3_cell = CellBasic(name='region 3')
    region3_cell.setMaterial(materials['region_3'])
    region3 = Universe(name='region 3')
    region3.addCell(region3_cell)

    # Region 4
    region4_cell = CellBasic(name='region 4')
    region4_cell.setMaterial(materials['region_4'])
    region4 = Universe(name='region 4')
    region4.addCell(region4_cell)

    # Region 5
    region5_cell = CellBasic(name='region 5')
    region5_cell.setMaterial(materials['region_5'])
    region5 = Universe(name='region 5')
    region5.addCell(region5_cell)

    # Region 5
    region6_cell = CellBasic(name='region 6')
    region6_cell.setMaterial(materials['region_6'])
    region6 = Universe(name='region 6')
    region6.addCell(region6_cell)

    # CellFills
    assembly1_cell = CellFill(name='assembly 1')
    assembly2_cell = CellFill(name='assembly 2')
    assembly3_cell = CellFill(name='assembly 3')
    assembly4_cell = CellFill(name='assembly 4')
    assembly5_cell = CellFill(name='assembly 5')
    assembly6_cell = CellFill(name='assembly 6')

    assembly1 = Universe(name='assembly 1')
    assembly2 = Universe(name='assembly 2')
    assembly3 = Universe(name='assembly 3')
    assembly4 = Universe(name='assembly 4')
    assembly5 = Universe(name='assembly 5')
    assembly6 = Universe(name='assembly 6')

    assembly1.addCell(assembly1_cell)
    assembly2.addCell(assembly2_cell)
    assembly3.addCell(assembly3_cell)
    assembly4.addCell(assembly4_cell)
    assembly5.addCell(assembly5_cell)
    assembly6.addCell(assembly6_cell)

    # Root cell/universe
    root_cell = CellFill(name='root cell')
    root_cell.addSurface(halfspace=+1, surface=left)
    root_cell.addSurface(halfspace=-1, surface=right)
    root_cell.addSurface(halfspace=+1, surface=bottom)
    root_cell.addSurface(halfspace=-1, surface=top)

    root_universe = Universe(name='root universe')
    root_universe.addCell(root_cell)


    ###############################################################################
    ###########################   Creating Lattices   #############################
    ###############################################################################

    log.py_printf('NORMAL', 'Creating LRA lattices...')

    # Assembly 1
    assembly1_lattice = Lattice(name='assembly 1')
    assembly1_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region1] * 10] * 10
    assembly1_lattice.setUniverses(template)
    assembly1_cell.setFill(assembly1_lattice)

    # Assembly 2
    assembly2_lattice = Lattice(name='assembly 2')
    assembly2_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region2] * 10] * 10
    assembly2_lattice.setUniverses(template)
    assembly2_cell.setFill(assembly2_lattice)

    # Assembly 3
    assembly3_lattice = Lattice(name='assembly 3')
    assembly3_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region3] * 10] * 10
    assembly3_lattice.setUniverses(template)
    assembly3_cell.setFill(assembly3_lattice)

    # Assembly 4
    assembly4_lattice = Lattice(name='assembly 4')
    assembly4_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region4] * 10] * 10
    assembly4_lattice.setUniverses(template)
    assembly4_cell.setFill(assembly4_lattice)

    # Assembly 5
    assembly5_lattice = Lattice(name='assembly 5')
    assembly5_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region5] * 10] * 10
    assembly5_lattice.setUniverses(template)
    assembly5_cell.setFill(assembly5_lattice)

    # Assembly 6
    assembly6_lattice = Lattice(name='assembly 6')
    assembly6_lattice.setWidth(width_x=1.5, width_y=1.5)
    template = [[region6] * 10] * 10
    assembly6_lattice.setUniverses(template)
    assembly6_cell.setFill(assembly6_lattice)

    # Full core
    core_lattice = Lattice(name='core')
    core_lattice.setWidth(width_x=15.0, width_y=15.0)

    universes = {7 : assembly1, 8 : assembly2, 9: assembly3,
                 10 : assembly4, 11 : assembly5, 12 : assembly6}
    template = [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                [ 9,  9,  9,  9,  9,  9,  9, 12, 12, 12, 12],
                [ 9,  9,  9,  9,  9,  9,  9, 10, 12, 12, 12],
                [ 8,  7,  7,  7,  7,  8,  8, 11, 11, 12, 12],
                [ 8,  7,  7,  7,  7,  8,  8, 11, 11, 12, 12],
                [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
                [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
                [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
                [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
                [ 8,  7,  7,  7,  7,  8,  8,  9,  9, 12, 12]]

    for i in range(11):
      for j in range(11):
        template[i][j] = universes[template[i][j]]
    core_lattice.setUniverses(template)
    root_cell.setFill(core_lattice)


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
    solver.setSourceConvergenceThreshold(tolerance)
    solver.setNumThreads(num_threads)
    solver.convergeSource(max_iters)
    Keff = solver.getKeff()
    return Keff

# assign values for use in test case instance
test_type = 'Keff'
benchmark = 'LRA'
benchmark_value = benchmark_val
print benchmark_value, 'is the Keff found for', benchmark
error_margin = 0.0001
filename = 'LRA.py'
#setup_func = setup_LRA

test_LRA_case = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, num_threads='DEFAULT')
test_LRA_1t = regression_test_case(test_type, benchmark, benchmark_value, error_margin, filename, num_threads=1)

test_list = [(test_LRA_case, __name__), (test_LRA_1t, __name__)]

if __name__ == '__main__':
    
    output = open("LRA_failed_results.txt", "a") # create output file in case of failures
    LRA_test_suite = regression_test_suite(test_list, output)
    LRA_test_suite.run_tests()

