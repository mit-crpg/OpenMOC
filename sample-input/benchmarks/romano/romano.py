import numpy
from openmoc import *
import openmoc.log as log
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

log.py_printf('TITLE', 'Simulating HW3 from Fall 2010 22.212...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

fuel = Material(name='fuel')
moderator = Material(name='moderator')

fuel.setNumEnergyGroups(1)
moderator.setNumEnergyGroups(1)

fuel.setSigmaA(numpy.array([0.069389522]))
fuel.setSigmaT(numpy.array([0.452648699]))
fuel.setSigmaF(numpy.array([0.0414198575]))
fuel.setNuSigmaF(numpy.array([0.0994076580]))
fuel.setSigmaS(numpy.array([0.38259177]))
fuel.setChi(numpy.array([1.0]))

moderator.setSigmaA(numpy.array([0.003751099]))
moderator.setSigmaT(numpy.array([0.841545641]))
moderator.setSigmaF(numpy.array([0.0]))
moderator.setNuSigmaF(numpy.array([0.0]))
moderator.setSigmaS(numpy.array([0.837794542]))
moderator.setChi(numpy.array([1.0]))


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=0.4)
left = XPlane(x=-0.635)
right = XPlane(x=0.635)
top = YPlane(y=0.635)
bottom = YPlane(y=-0.635)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel_cell = CellBasic(name='fuel')
fuel_cell.setMaterial(fuel)
fuel_cell.addSurface(halfspace=-1, surface=circle)

moderator_cell = CellBasic(name='moderator')
moderator_cell.setMaterial(moderator)
moderator_cell.addSurface(halfspace=+1, surface=circle)
moderator_cell.addSurface(halfspace=+1, surface=left)
moderator_cell.addSurface(halfspace=-1, surface=right)
moderator_cell.addSurface(halfspace=+1, surface=bottom)
moderator_cell.addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(fuel_cell)
root_universe.addCell(moderator_cell)


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

log.py_printf('TITLE', 'Finished')
