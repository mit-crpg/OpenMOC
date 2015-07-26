import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
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
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

sigma_a = numpy.array([0.069389522])
sigma_f = numpy.array([0.0414198575])
nu_sigma_f = numpy.array([0.0994076580])
sigma_s = numpy.array([0.383259177])
chi = numpy.array([1.0])
sigma_t = numpy.array([0.452648699])

infinite_medium = Material(name='1-group infinite medium')
infinite_medium.setNumEnergyGroups(1)
infinite_medium.setSigmaA(sigma_a)
infinite_medium.setSigmaF(sigma_f)
infinite_medium.setNuSigmaF(nu_sigma_f)
infinite_medium.setSigmaS(sigma_s)
infinite_medium.setChi(chi)
infinite_medium.setSigmaT(sigma_t)


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
#                             Creating Universes
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

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

log.py_printf('NORMAL', 'Running MOC eigenvalue simulation...')

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()


###############################################################################
#                            Verify Running a Simulation
###############################################################################

log.py_printf('NORMAL', 'Verifying with NumPy eigenvalue...')

# Compute fission production matrix
fiss_mat = numpy.outer(chi, nu_sigma_f)

# Create forward operator
M = numpy.linalg.solve((numpy.diag(sigma_t) - sigma_s), fiss_mat)

# Solve eigenvalue problem with NumPy
k, phi = numpy.linalg.eig(M)

# Select the dominant eigenvalue
k = k[0]

log.py_printf('RESULT', 'Numpy eigenvalue: {0:.6f}'.format(k))

log.py_printf('TITLE', 'Finished')
