## Regression Testing
## Author: Luisa Kenausis
## kenausis@mit.edu

## Should be run from within OpenMOC/sample-input/

from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options
import unittest
import sys

def generic_test_setup(sysargs):

    ## To clean up code: this sets up for a test and takes in the command
    ## line arguments that should be tested.

    sys.argv = sysargs

    options = Options()    

    num_threads = options.getNumThreads()
    track_spacing = options.getTrackSpacing()
    num_azim = options.getNumAzimAngles()
    tolerance = options.getTolerance()
    max_iters = options.getMaxIterations()

    log.set_log_level('ERROR')

    materials = materialize.materialize('c5g7-materials.h5')
    
    uo2_id = materials['UO2'].getId()
    water_id = materials['Water'].getId()

    circle = Circle(x=0.0, y=0.0, radius=1.0)
    left = XPlane(x=-2.0)
    right = XPlane(x=2.0)
    top = YPlane(y=2.0)
    bottom = YPlane(y=-2.0)

    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)

    ## Creating Cells

    cells = []
    cells.append(CellBasic(universe=1, material=uo2_id))
    cells.append(CellBasic(universe=1, material=water_id))
    cells.append(CellFill(universe=0, universe_fill=2))

    cells[0].addSurface(halfspace=-1, surface=circle)
    cells[1].addSurface(halfspace=+1, surface=circle)
    cells[2].addSurface(halfspace=+1, surface=left)
    cells[2].addSurface(halfspace=-1, surface=right)
    cells[2].addSurface(halfspace=+1, surface=bottom)
    cells[2].addSurface(halfspace=-1, surface=top)

    ## Creating Lattices

    lattice = Lattice(id=2, width_x=4.0, width_y=4.0)
    lattice.setLatticeCells([[1]])

    ## Creating Geometry

    g = Geometry()
    for material in materials.values(): g.addMaterial(material)
    for cell in cells: g.addCell(cell)
    g.addLattice(lattice)

    g.initializeFlatSourceRegions()

    ## Creating the TrackGenerator

    track_generator = TrackGenerator(g, num_azim, track_spacing)
    track_generator.generateTracks()

    ## Running a Simulation

    solver = ThreadPrivateSolver(g, track_generator)
    solver.setNumThreads(num_threads)
    solver.setSourceConvergenceThreshold(tolerance)
    Keff = solver.getKeff()
    TotalTime =  solver.getTotalTime()
    NumPolarAngles = solver.getNumPolarAngles()
    PolarQuadratureType = solver.getPolarQuadratureType()
    NumFSRs = g.getNumFSRs()
    SourceConvergenceThreshold = solver.getSourceConvergenceThreshold()
    NumGroups = g.getNumEnergyGroups()
    NumTracks = track_generator.getNumTracks()
    NumSegments = track_generator.getNumSegments()

    return [Keff, TotalTime, NumPolarAngles, PolarQuadratureType, NumFSRs,
            SourceConvergenceThreshold, NumGroups, NumTracks, NumSegments]
    

class TestRegression(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        ## run through solver with default settings except for
        ## max-iters at 25
        
        vals = generic_test_setup(["pin-cell.py", '--max-iters', '50'])

        cls._Keff = vals[0]
        cls._TotalTime = vals[1]
        cls._NumPolarAngles = vals[2]
        cls._PolarQuadratureType = vals[3]
        cls._NumFSRs = vals[4]
        cls._SourceConvergenceThreshold = vals[5]
        cls._NumGroups = vals[6]
        cls._NumTracks = vals[7]
        cls._NumSegments = vals[8]

    def setUp(self):

        vals = generic_test_setup(["pin-cell.py", '--max-iters', '50'])

        self.__Keff = vals[0]
        self.__TotalTime = vals[1]
        self.__NumPolarAngles = vals[2]
        self.__PolarQuadratureType = vals[3]
        self.__NumFSRs = vals[4]
        self.__SourceConvergenceThreshold = vals[5]
        self.__NumGroups = vals[6]
        self.__NumTracks = vals[7]
        self.__NumSegments = vals[8]

    def testKeff(self):

        ## Create a simulation w/ exact same setup
        ## and compare its calculated Keff with setup

        self.assertEqual(self._Keff, self.__Keff)

        ## This sometimes works (when both are found to be around 1.65e-33)
        ## but sometimes Keff is found to be on the order of 1e+27...

        ## I've duplicated it several times here so that the variability is
        ## more likely to be seen running it.

    def testKeff2(self):

        self.testKeff()

    def testKeff3(self):

        self.testKeff()
        
    def testKeff4(self):

        self.testKeff()
        
    def testKeff5(self):

        self.testKeff()
        
    def testTime(self):

        self.assertEqual(self._TotalTime, self.__TotalTime)

    def testNumPolarAngles(self):

        self.assertEqual(self._NumPolarAngles, self.__NumPolarAngles)

    def testPolarQuadratureType(self):

        self.assertEqual(self._PolarQuadratureType, self.__PolarQuadratureType)

    def testNumFSRs(self):

        self.assertEqual(self._NumFSRs, self.__NumFSRs)

    def testSourceConvergenceThreshold(self):

        self.assertEqual(self._SourceConvergenceThreshold, self.__SourceConvergenceThreshold)

    def testNumGroups(self):

        self.assertEqual(self._NumGroups, self.__NumGroups)

    def testNumTracks(self):

        self.assertEqual(self._NumTracks, self.__NumTracks)

    def testNumSegments(self):

        self.assertEqual(self._NumSegments, self.__NumSegments)

    @classmethod
    def tearDownClass(cls):

        sys.argv = ['pin-cell.py']
        
class TestVariedThreads(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        
        ## run simulation with default settings

        vals = generic_test_setup(["pin-cell.py", "--num-omp-threads", "3"])

        cls._Keff = vals[0]
        cls._TotalTime = vals[1]
        cls._NumPolarAngles = vals[2]
        cls._PolarQuadratureType = vals[3]
        cls._NumFSRs = vals[4]
        cls._SourceConvergenceThreshold = vals[5]
        cls._NumGroups = vals[6]
        cls._NumTracks = vals[7]
        cls._NumSegments = vals[8]

    def testMoreThreadsKeff(self):

        ## run through solver with default settings except for
        ## num-omp-threads at 10 instead of 3; compare Keff

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._Keff, vals[0])

    def testFewerThreadsKeff(self):

        ## run through solver with default settings except for
        ## num-omp-threads at 1 instead of 3; compare Keff

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._Keff, vals[0])

    def testMoreThreadsTime(self):

        ## run through solver with default settings except for
        ## num-omp-threads at 10 instead of 3; compare total time

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._TotalTime, vals[1])


    def testFewerThreadsTime(self):

        ## run through solver with default settings except for
        ## num-omp-threads at 1 instead of 3; compare total time

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._TotalTime, vals[1])

    def testMoreThreadsPolarAngles(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._NumPolarAngles, vals[2])

    def testFewerThreadsPolarAngles(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._NumPolarAngles, vals[2])

    def testMoreThreadsPolarQuadType(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._PolarQuadratureType, vals[3])

    def testFewerThreadsPolarQuadType(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._PolarQuadratureType, vals[3])

    def testMoreThreadsNumFSR(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._NumFSRs, vals[4])

    def testFewerThreadsNumFSR(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._NumFSRs, vals[4])

    def testMoreThreadsSCT(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._SourceConvergenceThreshold, vals[5])

    def testFewerThreadsSCT(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._SourceConvergenceThreshold, vals[5])

    def testMoreThreadsNumGroups(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._NumGroups, vals[6])

    def testFewerThreadsNumGroups(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._NumGroups, vals[6])

    def testMoreThreadsNumTracks(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._NumTracks, vals[7])

    def testFewerThreadsNumTracks(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._NumTracks, vals[7])

    def testMoreThreadsNumSegments(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '10'])

        self.assertEqual(self._NumSegments, vals[8])

    def testFewerThreadsNumSegments(self):

        vals = generic_test_setup(["pin-cell.py", '--num-omp-threads', '1'])

        self.assertEqual(self._NumSegments, vals[8])


## for test-writing ref:

##        self.__Keff = vals[0]
##        self.__TotalTime = vals[1]
##        self.__NumPolarAngles = vals[2]
##        self.__PolarQuadratureType = vals[3]
##        self.__NumFSRs = vals[4]
##        self.__SourceConvergenceThreshold = vals[5]
##        self.__NumGroups = vals[6]
##        self.__NumTracks = vals[7]
##        self.__NumSegments = vals[8]
##
        
suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestVariedThreads))

unittest.TextTestRunner(verbosity=2).run(suite)
