import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.cuda as cuda


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

num_blocks = 64
num_threads = 64
track_spacing = 0.1
num_azim = 16
tolerance = 1E-5
max_iters = 1000

log.setLogLevel('NORMAL')

log.py_printf('TITLE', 'Simulating a one group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.43...')


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

infinite_medium = Material(1)
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

circle = Circle(x=0.0, y=0.0, radius=10.0)
left = XPlane(x=-100.0)
right = XPlane(x=100.0)
top = YPlane(y=100.0)
bottom = YPlane(y=-100.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=1, rings=2, sectors=4))
cells.append(CellBasic(universe=1, material=1))
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

lattice = Lattice(id=2, width_x=200.0, width_y=200.0)
lattice.setLatticeCells([[1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

Timer.startTimer()

geometry = Geometry()
geometry.addMaterial(infinite_medium)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()

Timer.stopTimer()
Timer.recordSplit('Iniitilializing the geometry')
Timer.resetTimer()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

Timer.startTimer()

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()

Timer.stopTimer()
Timer.recordSplit('Ray tracing across the geometry')
Timer.resetTimer()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

Timer.startTimer()

solver = cuda.GPUSolver(geometry, track_generator)
solver.setNumThreadBlocks(num_blocks)
solver.setNumThreadsPerBlock(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)

Timer.stopTimer()
Timer.recordSplit('Converging the source on the GPU')
Timer.resetTimer()
Timer.printSplits()

log.py_printf('TITLE', 'Finished')
