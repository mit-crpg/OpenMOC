import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('INFO')

log.py_printf('TITLE', 'Simulating a one group homogeneous infinite medium...')
log.py_printf('HEADER', 'The reference keff = 1.43...')

log.py_printf('NORMAL', 'Creating materials...need to use hdf5...')

infinite_medium = Material(1)
infinite_medium.setNumEnergyGroups(1)
infinite_medium.setSigmaA(numpy.array([0.069389522]))
infinite_medium.setSigmaF(numpy.array([0.0414198575]))
infinite_medium.setNuSigmaF(numpy.array([0.0994076580]))
infinite_medium.setSigmaS(numpy.array([0.383259177]))
infinite_medium.setChi(numpy.array([1.0]))
infinite_medium.setSigmaT(numpy.array([0.452648699]))

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(id=1, x=0.0, y=0.0, radius=10.0)
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

geometry = Geometry()
geometry.addMaterial(infinite_medium)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator()
track_generator.setNumAzim(128)
track_generator.setTrackSpacing(0.5)
track_generator.setGeometry(geometry)
track_generator.generateTracks()



solver = Solver(geometry, track_generator)
solver.setNumThreads(1)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(1000)
