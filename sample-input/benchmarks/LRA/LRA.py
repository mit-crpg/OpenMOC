from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

tolerance = 1E-10
log.set_log_level('DEBUG')

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data...')

materials = materialize.materialize('LRA-materials.py')

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
cells.append(CellFill(universe=0, universe_fill=7))

cells[6].addSurface(halfspace=+1, surface=planes[0])
cells[6].addSurface(halfspace=-1, surface=planes[1])
cells[6].addSurface(halfspace=+1, surface=planes[2])
cells[6].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattice...')

lattice = Lattice(id=7, width_x=15.0, width_y=15.0)
lattice.setLatticeCells([[6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                         [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                         [3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6],
                         [3, 3, 3, 3, 3, 3, 3, 4, 6, 6, 6],
                         [2, 1, 1, 1, 1, 2, 2, 5, 5, 6, 6],
                         [2, 1, 1, 1, 1, 2, 2, 5, 5, 6, 6],
                         [1, 1, 1, 1, 1, 1, 1, 3, 3, 6, 6],
                         [1, 1, 1, 1, 1, 1, 1, 3, 3, 6, 6],
                         [1, 1, 1, 1, 1, 1, 1, 3, 3, 6, 6],
                         [1, 1, 1, 1, 1, 1, 1, 3, 3, 6, 6],
                         [2, 1, 1, 1, 1, 2, 2, 3, 3, 6, 6]])


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
geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd = Cmfd(geometry)
cmfd.computeKeff()

log.py_printf('NORMAL', 'k_eff = %f', cmfd.getKeff())

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
plotter.plot_mesh_fluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

