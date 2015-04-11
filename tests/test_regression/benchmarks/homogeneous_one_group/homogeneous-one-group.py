import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
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

log.py_printf('TITLE', 'Simulating a one group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.43...')


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

infinite_medium = Material(name='1-group infinite medium')
infinite_medium.setNumEnergyGroups(1)
infinite_medium.setSigmaA(numpy.array([0.069389522]))
infinite_medium.setSigmaF(numpy.array([0.0414198575]))
infinite_medium.setNuSigmaF(numpy.array([0.0994076580]))
infinite_medium.setSigmaS(numpy.array([0.383259177]))
infinite_medium.setChi(numpy.array([1.0]))
infinite_medium.setSigmaT(numpy.array([0.452648699]))


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-100.0, name='left')
right = XPlane(x=100.0, name='right')
top = YPlane(y=100.0, name='top')
bottom = YPlane(y=-100.0, name='bottom')

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cell = CellBasic()
cell.setMaterial(infinite_medium)
cell.addSurface(halfspace=+1, surface=left)
cell.addSurface(halfspace=-1, surface=right)
cell.addSurface(halfspace=+1, surface=bottom)
cell.addSurface(halfspace=-1, surface=top)


###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(cell)


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
