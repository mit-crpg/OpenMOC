## test homogeneous-one-group

## reliably works + passes

import numpy
from openmoc import *
import openmoc.log as log
from openmoc.options import Options
import unittest
import sys

def general_h2g_setup(sysargs):
    
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

class test_h2g(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls._Keff_benchmark = 1.723158836364746
        cls._Keff = general_h2g_setup(['homogeneous-two-groups.py'])

    def test_h2g_Keff(self):

        self.assertEqual(self._Keff, self._Keff_benchmark)

    def test_h2g_Keff_more_threads(self):

        Keff = general_h2g_setup(['homogeneous-one-group.py', '--num-omp-threads', '5'])
        self.assertEqual(Keff, self._Keff_benchmark)

        ## this fails, although it's close.
        ## use small threshold instead?

test_h2g = unittest.TestLoader().loadTestsFromTestCase(test_h2g)

if __name__ == '__main__':    
    unittest.TextTestRunner(verbosity=2).run(test_h2g)

