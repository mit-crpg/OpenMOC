import numpy
from openmoc import *
import openmoc.log as log


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

num_threads = options.num_omp_threads
track_spacing = options.track_spacing
num_azim = options.num_azim
tolerance = options.tolerance
max_iters = options.max_iters

log.setLogLevel('NORMAL')

log.py_printf('TITLE', 'Simulating HW3 from Fall 2010 22.212...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

fuel = Material(1)
moderator = Material(2)

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

cells = []
cells.append(CellBasic(universe=1, material=1))
cells.append(CellBasic(universe=1, material=2))
cells.append(CellFill(universe=0, universe_fill=2))

cells[0].addSurface(halfspace=-1, surface=circle)
cells[1].addSurface(halfspace=+1, surface=circle)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=2, width_x=1.27, width_y=1.27)
lattice.setLatticeCells([[1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.addMaterial(fuel)
geometry.addMaterial(moderator)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = ThreadPrivateSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()
