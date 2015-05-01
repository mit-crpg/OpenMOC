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

circle = Circle(x=0.0, y=0.0, radius=1.0, name='pin')
left = XPlane(x=-2.0, name='left')
right = XPlane(x=2.0, name='right')
top = YPlane(y=2.0, name='top')
bottom = YPlane(y=-2.0, name='bottom')

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = CellBasic(name='fuel')
fuel.setMaterial(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=circle)

moderator = CellBasic(name='moderator')
moderator.setMaterial(materials['Water'])
moderator.addSurface(halfspace=+1, surface=circle)
moderator.addSurface(halfspace=+1, surface=left)
moderator.addSurface(halfspace=-1, surface=right)
moderator.addSurface(halfspace=+1, surface=bottom)
moderator.addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(fuel)
root_universe.addCell(moderator)


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

solver = CPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_segments(track_generator)
#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
#plotter.plot_energy_fluxes(solver, fsrs=range(geometry.getNumFSRs()))

log.py_printf('TITLE', 'Finished')
