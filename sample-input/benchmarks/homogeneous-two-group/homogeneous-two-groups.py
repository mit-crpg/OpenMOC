from openmoc import *
import numpy
import openmoc.log as log
import openmoc.plotter as plotter


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

num_threads = 4
track_spacing = 0.1
num_azim = 16
tolerance = 1E-3
max_iters = 1000

log.py_setlevel('INFO')

log.py_printf('TITLE', 'Simulating a two group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.72...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

infinite_medium = Material(1)
infinite_medium.setNumEnergyGroups(2)
infinite_medium.setSigmaA(numpy.array([0.0038, 0.184]))
infinite_medium.setSigmaF(numpy.array([0.000625, 0.135416667]))
infinite_medium.setNuSigmaF(numpy.array([0.0015, 0.325]))
infinite_medium.setSigmaS(numpy.array([0.1, 0.117, 0.0, 1.42]))
infinite_medium.setChi(numpy.array([1.0, 0.0]))
infinite_medium.setSigmaT(numpy.array([0.2208, 1.604]))


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=50.0)
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
cells.append(CellBasic(universe=1, material=1))
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

geometry = Geometry()
geometry.addMaterial(infinite_medium)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator()
track_generator.setNumAzim(num_azim)
track_generator.setTrackSpacing(track_spacing)
track_generator.setGeometry(geometry)
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = Solver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)

Timer.startTimer()
solver.convergeSource(max_iters)
Timer.stopTimer()
Timer.recordSplit('Fixed source iteration on host')
Timer.printSplits()

log.py_printf('TITLE', 'Finished')
