import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('INFO')

log.py_printf('TITLE', 'Simulating HW3 from Fall 2010 22.212...')

log.py_printf('NORMAL', 'Creating materials...need to use hdf5...')

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

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(id=1, x=0.0, y=0.0, radius=0.4)
left = XPlane(id=2, x=-0.635)
right = XPlane(id=3, x=0.635)
top = YPlane(id=4, y=0.635)
bottom = YPlane(id=5, y=-0.635)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(id=1, universe=1, material=1))
cells.append(CellBasic(id=2, universe=1, material=2))
cells.append(CellFill(id=3, universe=0, universe_fill=2))

cells[0].addSurface(halfspace=-1, surface=circle)
cells[1].addSurface(halfspace=+1, surface=circle)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=2, width_x=1.27, width_y=1.27)
lattice.setLatticeCells([[1]])

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.addMaterial(fuel)
geometry.addMaterial(moderator)
geometry.addCell(cells[0])
geometry.addCell(cells[1])
geometry.addCell(cells[2])
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator()
track_generator.setNumAzim(64)
track_generator.setTrackSpacing(0.1)
track_generator.setGeometry(geometry)
track_generator.generateTracks()

solver = Solver(geometry, track_generator)
solver.setNumThreads(1)
solver.setSourceConvergenceThreshold(1E-4)
solver.convergeSource(10000)
