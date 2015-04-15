from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-10.0, name='left')
right = XPlane(x=10.0, name='right')
top = YPlane(y=-10.0, name='top')
bottom = YPlane(y=10.0, name='bottom')
boundaries = [left, right, top, bottom]

for boundary in boundaries: boundary.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

moderator_cell = CellBasic(name='moderator')
moderator_cell.setMaterial(Material())

root_cell = CellFill(name='root cell')
root_cell.addSurface(halfspace=+1, surface=boundaries[0])
root_cell.addSurface(halfspace=-1, surface=boundaries[1])
root_cell.addSurface(halfspace=+1, surface=boundaries[2])
root_cell.addSurface(halfspace=-1, surface=boundaries[3])


###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

moderator_univ = Universe(name='moderator')
root_universe = Universe(name='root universe')

moderator_univ.addCell(moderator_cell)
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 4 x 4 lattice...')

num_x = 10
num_y = 10
width_x = (root_universe.getMaxX() - root_universe.getMinX()) / num_y
width_y = (root_universe.getMaxY() - root_universe.getMinY()) / num_x
universes = [[moderator_univ]*num_x]*num_y

lattice = Lattice(name='{0}x{1} lattice'.format(num_x, num_y))
lattice.setWidth(width_x=width_x, width_y=width_y)
lattice.setUniverses(universes)
root_cell.setFill(lattice)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)

solver.setFixedSourceByFSR(25, 7, 1.0)

solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=250)
plotter.plot_cells(geometry, gridsize=250)
plotter.plot_flat_source_regions(geometry, gridsize=250)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
