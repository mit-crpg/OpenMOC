from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.set_log_level('INFO')

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from py...')

materials = materialize.materialize('materialsLRA.py')

region1 = materials['region_1'].getId()
region2 = materials['region_2'].getId()
region3 = materials['region_3'].getId()
region4 = materials['region_4'].getId()
region5 = materials['region_5'].getId()
region6 = materials['region_6'].getId()

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

planes = []
planes.append(XPlane(x=-82.5))
planes.append(XPlane(x=82.5))
planes.append(YPlane(y=-82.5))
planes.append(YPlane(y=82.5))
planes[0].setBoundaryType(REFLECTIVE)
planes[1].setBoundaryType(ZERO_FLUX)
planes[2].setBoundaryType(REFLECTIVE)
planes[3].setBoundaryType(ZERO_FLUX)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

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


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattice...')

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


###############################################################################
###########################   Creating Cmfd Mesh   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd mesh...')

mesh = Mesh(DIFFUSION)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

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

geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd = Cmfd(geometry)
cmfd.setOmega(1.5)
cmfd.computeKeff()

log.py_printf('NORMAL', 'k_eff = %f', cmfd.getKeff())

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_mesh_fluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

