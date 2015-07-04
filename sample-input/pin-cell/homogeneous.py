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
azim_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()
num_polar = 6
polar_spacing = 0.1

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

xmin = XPlane(x=-0.5, name='xmin')
xmax = XPlane(x= 0.5, name='xmax')
ymin = YPlane(y=-0.5, name='ymin')
ymax = YPlane(y= 0.5, name='ymax')
zmin = ZPlane(z=-0.5, name='zmin')
zmax = ZPlane(z= 0.5, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

moderator = Cell(name='moderator')
moderator.setFill(materials['UO2'])
moderator.addSurface(halfspace=+1, surface=xmin)
moderator.addSurface(halfspace=-1, surface=xmax)
moderator.addSurface(halfspace=+1, surface=ymin)
moderator.addSurface(halfspace=-1, surface=ymax)
moderator.addSurface(halfspace=+1, surface=zmin)
moderator.addSurface(halfspace=-1, surface=zmax)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
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

quad = EqualAnglePolarQuad()
quad.setNumPolarAngles(num_polar)

track_generator = TrackGenerator(geometry, num_azim, num_polar, azim_spacing, \
                                 polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
#track_generator.setSolve2D()
track_generator.setTrackGenerationMethod(MODULAR_RAY_TRACING)
track_generator.generateTracks()

###############################################################################
#                            Running a Simulation
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

#plotter.plot_tracks(track_generator)
#plotter.plot_tracks_3d(track_generator)
#plotter.plot_segments_3d(track_generator)
#plotter.plot_segments(track_generator)
#plotter.plot_quadrature(track_generator.getQuadrature())

log.py_printf('TITLE', 'Finished')
