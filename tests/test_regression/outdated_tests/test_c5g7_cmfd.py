## Tests c5g7-cmfd benchmark for correct K_eff value.

## Sample Test for Code Demo -- not as developed as most of the
##  existing benchmark tests.

from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options
import unittest
import sys
import os
import time

## This adds tests/test_regression/ to the path list so it can access and import
## the run_regression_test module.

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_directory[:-21])

from regression_test_runner import run_regression_test, write_failed_results

testfail = False # variable to use when creating test file
output = open("c5g7_cmfd_failed_results2.txt", "a+") # create output file in case of failures

def general_c5g7_cmfd_setup(sysargs):

    sys.argv = sysargs

    ## parameters
    options = Options()

    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()
    log.set_log_level('ERROR')

    ## materials
    materials = materialize.materialize(current_directory + '/c5g7_materials.h5')


## To Break -- change materials['UO2']
    ## setNuSigmaFByGroup(val, group)        

#    materials['UO2'].setNuSigmaFByGroup(materials['UO2'].getSigmaFByGroup(2)*.9,2)
    uo2_id = materials['UO2'].getId()
    mox43_id = materials['MOX-4.3%'].getId()
    mox7_id = materials['MOX-7%'].getId()
    mox87_id = materials['MOX-8.7%'].getId()
    guide_tube_id = materials['Guide Tube'].getId()
    fiss_id = materials['Fission Chamber'].getId()
    water_id = materials['Water'].getId()
        
    ## surfaces

    circles = []
    planes = []
    planes.append(XPlane(x=-32.13))
    planes.append(XPlane(x=32.13))
    planes.append(YPlane(y=-32.13))
    planes.append(YPlane(y=32.13))
    circles.append(Circle(x=0., y=0., radius=0.54))
    circles.append(Circle(x=0., y=0., radius=0.58))
    circles.append(Circle(x=0., y=0., radius=0.62))
    planes[0].setBoundaryType(REFLECTIVE)
    planes[1].setBoundaryType(VACUUM)
    planes[2].setBoundaryType(VACUUM)
    planes[3].setBoundaryType(REFLECTIVE)

    ## cells
    cells = []

    # UO2 pin cells
    cells.append(CellBasic(universe=1, material=uo2_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=1, material=water_id, sectors=8))
    cells.append(CellBasic(universe=1, material=water_id, sectors=8))
    cells.append(CellBasic(universe=1, material=water_id, sectors=8))
    cells[0].addSurface(-1, circles[0])
    cells[1].addSurface(+1, circles[0])
    cells[1].addSurface(-1, circles[1])
    cells[2].addSurface(+1, circles[1])
    cells[2].addSurface(-1, circles[2])
    cells[3].addSurface(+1, circles[2])


    # 4.3% MOX pin cells
    cells.append(CellBasic(universe=2, material=mox43_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=2, material=water_id, sectors=8))
    cells.append(CellBasic(universe=2, material=water_id, sectors=8))
    cells.append(CellBasic(universe=2, material=water_id, sectors=8))
    cells[4].addSurface(-1, circles[0])
    cells[5].addSurface(+1, circles[0])
    cells[5].addSurface(-1, circles[1])
    cells[6].addSurface(+1, circles[1])
    cells[6].addSurface(-1, circles[2])
    cells[7].addSurface(+1, circles[2])


    # 7% MOX pin cells
    cells.append(CellBasic(universe=3, material=mox7_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=3, material=water_id, sectors=8))
    cells.append(CellBasic(universe=3, material=water_id, sectors=8))
    cells.append(CellBasic(universe=3, material=water_id, sectors=8))
    cells[8].addSurface(-1, circles[0])
    cells[9].addSurface(+1, circles[0])
    cells[9].addSurface(-1, circles[1])
    cells[10].addSurface(+1, circles[1])
    cells[10].addSurface(-1, circles[2])
    cells[11].addSurface(+1, circles[2])


    # 8.7% MOX pin cells
    cells.append(CellBasic(universe=4, material=mox87_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=4, material=water_id, sectors=8))
    cells.append(CellBasic(universe=4, material=water_id, sectors=8))
    cells.append(CellBasic(universe=4, material=water_id, sectors=8))
    cells[12].addSurface(-1, circles[0])
    cells[13].addSurface(+1, circles[0])
    cells[13].addSurface(-1, circles[1])
    cells[14].addSurface(+1, circles[1])
    cells[14].addSurface(-1, circles[2])
    cells[15].addSurface(+1, circles[2])

    # Fission chamber pin cells
    cells.append(CellBasic(universe=5, material=fiss_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=5, material=water_id, sectors=8))
    cells.append(CellBasic(universe=5, material=water_id, sectors=8))
    cells.append(CellBasic(universe=5, material=water_id, sectors=8))
    cells[16].addSurface(-1, circles[0])
    cells[17].addSurface(+1, circles[0])
    cells[17].addSurface(-1, circles[1])
    cells[18].addSurface(+1, circles[1])
    cells[18].addSurface(-1, circles[2])
    cells[19].addSurface(+1, circles[2])

    # Guide tube pin cells
    cells.append(CellBasic(universe=6, material=guide_tube_id, rings=3, sectors=8))
    cells.append(CellBasic(universe=6, material=water_id, sectors=8))
    cells.append(CellBasic(universe=6, material=water_id, sectors=8))
    cells.append(CellBasic(universe=6, material=water_id, sectors=8))
    cells[20].addSurface(-1, circles[0])
    cells[21].addSurface(+1, circles[0])
    cells[21].addSurface(-1, circles[1])
    cells[22].addSurface(+1, circles[1])
    cells[22].addSurface(-1, circles[2])
    cells[23].addSurface(+1, circles[2])

    # Moderator cell
    cells.append(CellBasic(universe=7, material=water_id))

    # Top left, bottom right lattice
    cells.append(CellFill(universe=10, universe_fill=20))

    # Top right, bottom left lattice
    cells.append(CellFill(universe=11, universe_fill=21))

    # Moderator lattice - semi-finely spaced
    cells.append(CellFill(universe=12, universe_fill=23))

    # Moderator lattice - bottom of geometry
    cells.append(CellFill(universe=13, universe_fill=24))

    # Moderator lattice - bottom corner of geometry
    cells.append(CellFill(universe=14, universe_fill=25))

    # Moderator lattice right side of geometry
    cells.append(CellFill(universe=15, universe_fill=26))

    # Full geometry
    cells.append(CellFill(universe=0, universe_fill=30))
    cells[-1].addSurface(+1, planes[0])
    cells[-1].addSurface(-1, planes[1])
    cells[-1].addSurface(+1, planes[2])
    cells[-1].addSurface(-1, planes[3])

    ## lattices
    lattices = []

    # Top left, bottom right 17 x 17 assemblies
    lattices.append(Lattice(id=20, width_x=1.26, width_y=1.26))
    lattices[-1].setLatticeCells(
        [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
         [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
         [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])


    # Top right, bottom left 17 x 17 assemblies 
    lattices.append(Lattice(id=21, width_x=1.26, width_y=1.26))
    lattices[-1].setLatticeCells(
        [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
         [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
         [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
         [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
         [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
         [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
         [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
         [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
         [2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2],
         [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
         [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
         [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
         [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
         [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
         [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
         [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
         [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]])


    # Sliced up water cells - semi finely spaced
    lattices.append(Lattice(id=23, width_x=0.126, width_y=0.126))
    lattices[-1].setLatticeCells(
        [[7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])


    # Sliced up water cells - right side of geometry
    lattices.append(Lattice(id=26, width_x=1.26, width_y=1.26))
    lattices[-1].setLatticeCells(
        [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7]])


    # Sliced up water cells for bottom corner of geometry
    lattices.append(Lattice(id=25, width_x=1.26, width_y=1.26))
    lattices[-1].setLatticeCells(
        [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])


    # Sliced up water cells for bottom of geometry
    lattices.append(Lattice(id=24, width_x=1.26, width_y=1.26))
    lattices[-1].setLatticeCells(
        [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
         [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])


    # 4 x 4 core to represent two bundles and water
    lattices.append(Lattice(id=30, width_x=21.42, width_y=21.42))
    lattices[-1].setLatticeCells([[10, 11, 15],
                                 [11, 10, 15],


                                 [13, 13, 14]])

    ## Cmfd mesh
    cmfd = Cmfd()
    cmfd.setMOCRelaxationFactor(0.6)
    cmfd.setSORRelaxationFactor(1.5)
    cmfd.setLatticeStructure(51,51)
    cmfd.setGroupStructure([1,4,8])

    ## geometry
    geometry = Geometry()
    geometry.setCmfd(cmfd)
    for material in materials.values(): geometry.addMaterial(material)
    for cell in cells: geometry.addCell(cell)
    for lattice in lattices: geometry.addLattice(lattice)

    geometry.initializeFlatSourceRegions()

    ## TrackGenerator
    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.setNumThreads(num_threads)
    track_generator.generateTracks()

    ## run simulation
    solver = CPUSolver(geometry, track_generator)
    solver.setSourceConvergenceThreshold(tolerance)
    solver.setNumThreads(num_threads)
    solver.convergeSource(max_iters)
    Keff = solver.getKeff()

    return Keff

class test_c5g7_cmfd(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # known benchmark value
        cls._Keff_benchmark = 1.1853487491607666

        ## run simulation
        start = time.clock()
        cls._Keff = general_c5g7_cmfd_setup(['c5g7-cmfd.py'])
        print 'finding Keff took ', time.clock()-start
##        self._test_type = 'Keff'
##        self._benchmark = 'c5g7_cmfd'
##        self._error_margin = 0.00000


    def test_Keff(self):
        self.assertAlmostEqual(self._Keff, self._Keff_benchmark, 3)

    def test_Keff_1t(self):
        Keff_1t = general_c5g7_cmfd_setup(['c5g7-cmfd.py', '--num-omp-threads', '1'])
        self.assertAlmostEqual(Keff_1t, self._Keff_benchmark, 3)

test_c5g7_cmfd = unittest.TestLoader().loadTestsFromTestCase(test_c5g7_cmfd)

if __name__ == '__main__':    
    unittest.TextTestRunner(verbosity=3).run(test_c5g7_cmfd)
    if testfail:
        output.write("------------------------------------------------------"+"\n")
    output.close()

