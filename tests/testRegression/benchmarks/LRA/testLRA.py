## Test LRA benchmark

## This gets different values sometimes (not reliable)

## ALSO: cmfd.computeKeff() has error

from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
import unittest
import sys

def general_LRA_setup(sysargs):

    ## Run simulation with passed-in parameters.

    sys.argv = sysargs
    log.set_log_level('ERROR')

    ## materials
    sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/LRA')
    materials = materialize.materialize('LRA-materials.py')

    ## CURRENTLY THROWS ERROR IN FULL TEST SUITE - 'no such file' when materializing
    ## 'LRA-materials.py' (could be hyphen?)

    region1 = materials['region_1'].getId()
    region2 = materials['region_2'].getId()
    region3 = materials['region_3'].getId()
    region4 = materials['region_4'].getId()
    region5 = materials['region_5'].getId()
    region6 = materials['region_6'].getId()

    ## surfaces
    planes = []
    planes.append(XPlane(x=-82.5))
    planes.append(XPlane(x=82.5))
    planes.append(YPlane(y=-82.5))
    planes.append(YPlane(y=82.5))
    planes[0].setBoundaryType(REFLECTIVE)
    planes[1].setBoundaryType(ZERO_FLUX)
    planes[2].setBoundaryType(REFLECTIVE)
    planes[3].setBoundaryType(ZERO_FLUX)

    ## cells
    cells = []
    cells.append(CellBasic(universe=1, material=region1))
    cells.append(CellBasic(universe=2, material=region2))
    cells.append(CellBasic(universe=3, material=region3))
    cells.append(CellBasic(universe=4, material=region4))
    cells.append(CellBasic(universe=5, material=region5))
    cells.append(CellBasic(universe=6, material=region6))
    cells.append(CellFill(universe=21, universe_fill=31))
    cells.append(CellFill(universe=22, universe_fill=32))
    cells.append(CellFill(universe=23, universe_fill=33))
    cells.append(CellFill(universe=24, universe_fill=34))
    cells.append(CellFill(universe=25, universe_fill=35))
    cells.append(CellFill(universe=26, universe_fill=36))
    cells.append(CellFill(universe=0, universe_fill=7))

    cells[12].addSurface(halfspace=+1, surface=planes[0])
    cells[12].addSurface(halfspace=-1, surface=planes[1])
    cells[12].addSurface(halfspace=+1, surface=planes[2])
    cells[12].addSurface(halfspace=-1, surface=planes[3])

    ## lattices
    assembly1 = Lattice(id=31, width_x=1.5, width_y=1.5)
    assembly1.setLatticeCells([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])

    assembly2 = Lattice(id=32, width_x=1.5, width_y=1.5)
    assembly2.setLatticeCells([[2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                               [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]])

    assembly3 = Lattice(id=33, width_x=1.5, width_y=1.5)
    assembly3.setLatticeCells([[3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                               [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]])


    assembly4 = Lattice(id=34, width_x=1.5, width_y=1.5)
    assembly4.setLatticeCells([[4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
                               [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]])

    assembly5 = Lattice(id=35, width_x=1.5, width_y=1.5)
    assembly5.setLatticeCells([[5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]])


    assembly6 = Lattice(id=36, width_x=1.5, width_y=1.5)
    assembly6.setLatticeCells([[6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                               [6, 6, 6, 6, 6, 6, 6, 6, 6, 6]])


    core = Lattice(id=7, width_x=15.0, width_y=15.0)
    core.setLatticeCells([[26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26],
                             [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26],
                             [23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26],
                             [23, 23, 23, 23, 23, 23, 23, 24, 26, 26, 26],
                             [22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26],
                             [22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26],
                             [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                             [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                             [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                             [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                             [22, 21, 21, 21, 21, 22, 22, 23, 23, 26, 26]])

    ## mesh
    mesh = Mesh(DIFFUSION)

    ## geometry
    geometry = Geometry(mesh)
    for material in materials.values(): geometry.addMaterial(material)
    for cell in cells: geometry.addCell(cell)
    geometry.addLattice(assembly1)
    geometry.addLattice(assembly2)
    geometry.addLattice(assembly3)
    geometry.addLattice(assembly4)
    geometry.addLattice(assembly5)
    geometry.addLattice(assembly6)
    geometry.addLattice(core)

    ## Cmfd module
    cmfd = Cmfd(geometry)
    cmfd.setOmega(1.5)
    ## cmfd.computeKeff()

    ## return Keff

    return cmfd.getKeff()

class TestLRA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls._Keff_benchmark = 1.04e-322
        cls._Keff = general_LRA_setup(['LRA.py'])

    def testKeff(self):

        self.assertEqual(self._Keff, self._Keff_benchmark)

    def testMoreThreadsKeff(self):

        Keff = general_LRA_setup(['LRA.py', '--num-omp-threads', '6'])
        self.assertEqual(Keff, self._Keff_benchmark)

suite = unittest.TestLoader().loadTestsFromTestCase(TestLRA)

unittest.TextTestRunner(verbosity=2).run(suite)

