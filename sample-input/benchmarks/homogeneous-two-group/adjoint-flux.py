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

log.py_printf('TITLE', 'Simulating a two group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.72...')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

sigma_a = numpy.array([0.0038, 0.184])
sigma_f = numpy.array([0.000625, 0.135416667])
nu_sigma_f = numpy.array([0.0015, 0.325])
sigma_s = numpy.array([[0.1, 0.117], [0., 1.42]])
chi = numpy.array([1.0, 0.0])
sigma_t = numpy.array([0.2208, 1.604])

infinite_medium = Material(name='2-group infinite medium')
infinite_medium.setNumEnergyGroups(2)
infinite_medium.setSigmaA(sigma_a)
infinite_medium.setSigmaF(sigma_f)
infinite_medium.setNuSigmaF(nu_sigma_f)
infinite_medium.setSigmaS(sigma_s.flat)
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

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

log.py_printf('NORMAL', 'Running MOC adjoint eigenvalue simulation...')

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters, mode=ADJOINT)
solver.printTimerReport()


###############################################################################
#                            Verify with NumPy
###############################################################################

log.py_printf('NORMAL', 'Verifying with NumPy adjoint eigenvalue...')

# Compute fission production matrix
fiss_mat = numpy.outer(chi, nu_sigma_f)

# Tranpose the scattering matrix from OpenMOC's notation
sigma_s = numpy.transpose(sigma_s)

# Create adjoint operator
sigma_s = numpy.transpose(sigma_s)
fiss_mat = numpy.transpose(fiss_mat)
M = numpy.linalg.solve((numpy.diag(sigma_t) - sigma_s), fiss_mat)

# Solve eigenvalue problem with NumPy
k, phi = numpy.linalg.eig(M)

# Select the dominant eigenvalue
k = max(k)

log.py_printf('RESULT', 'Numpy adjoint eigenvalue: {0:.6f}'.format(k))
log.py_printf('TITLE', 'Finished')
