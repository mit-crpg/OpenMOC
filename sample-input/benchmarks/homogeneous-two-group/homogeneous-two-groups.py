import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.num_omp_threads
azim_spacing = options.azim_spacing
num_azim = options.num_azim
polar_spacing = options.polar_spacing
num_polar = options.num_polar
tolerance = options.tolerance
max_iters = options.max_iters

log.set_log_level('NORMAL')

log.py_printf('TITLE', 'Simulating a two group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.72...')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

infinite_medium = Material(name='2-group infinite medium')
infinite_medium.setNumEnergyGroups(2)
#infinite_medium.setSigmaA(numpy.array([0.0038, 0.184]))
infinite_medium.setSigmaF(numpy.array([0.000625, 0.135416667]))
infinite_medium.setNuSigmaF(numpy.array([0.0015, 0.325]))
infinite_medium.setSigmaS(numpy.array([0.1, 0.117, 0.0, 1.42]))
infinite_medium.setChi(numpy.array([1.0, 0.0]))
infinite_medium.setSigmaT(numpy.array([0.2208, 1.604]))


###############################################################################
#                            Creating Surfaces
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
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cell = Cell()
cell.setFill(infinite_medium)
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
#                         Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, num_azim, azim_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

log.py_printf('TITLE', 'Finished')
