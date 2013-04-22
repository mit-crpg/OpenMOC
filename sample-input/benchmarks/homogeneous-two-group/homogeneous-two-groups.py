import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('INFO')

log.py_printf('TITLE', 'Simulating a two group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.72...')

timer = Timer()

log.py_printf('NORMAL', 'Creating materials...need to use hdf5...')

infinite_medium = Material(1)
infinite_medium.setNumEnergyGroups(2)
infinite_medium.setSigmaA(numpy.array([0.0038, 0.184]))
infinite_medium.setSigmaF(numpy.array([0.000625, 0.135416667]))
infinite_medium.setNuSigmaF(numpy.array([0.0015, 0.325]))
infinite_medium.setSigmaS(numpy.array([0.1, 0.117, 0.0, 1.42]))
infinite_medium.setChi(numpy.array([1.0, 0.0]))
infinite_medium.setSigmaT(numpy.array([0.2208, 1.604]))

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(id=1, x=0.0, y=0.0, radius=50.0)
left = XPlane(id=2, x=-100.0)
right = XPlane(id=3, x=100.0)
top = YPlane(id=4, y=100.0)
bottom = YPlane(id=5, y=-100.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(id=1, universe=1, material=1))
cells.append(CellBasic(id=2, universe=1, material=1))
cells.append(CellFill(id=3, universe=0, universe_fill=2))

cells[0].addSurface(halfspace=-1, surface=circle)
cells[1].addSurface(halfspace=+1, surface=circle)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=2, width_x=200.0, width_y=200.0)
lattice.setLatticeCells([[1]])

log.py_printf('NORMAL', 'Creating geometry...')

timer.startTimer()
geometry = Geometry()
geometry.addMaterial(infinite_medium)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()
timer.stopTimer()
timer.recordSplit('Geometry initialization')
timer.resetTimer()

log.py_printf('NORMAL', 'Initializing the track generator...')

timer.startTimer()
track_generator = TrackGenerator()
track_generator.setNumAzim(128)
track_generator.setTrackSpacing(0.1)
track_generator.setGeometry(geometry)
track_generator.generateTracks()
timer.stopTimer()
timer.recordSplit('Generating tracks')
timer.resetTimer()

timer.startTimer()
solver = Solver(geometry, track_generator)
solver.setNumThreads(8)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(1000)
timer.stopTimer()
timer.recordSplit('Fixed source iteration on host')
timer.printSplits()
