import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('INFO')

log.py_printf('NORMAL', 'Creating materials...need to use hdf5...')

fuel = Material(1)
moderator = Material(2)

fuel.setNumEnergyGroups(7)
moderator.setNumEnergyGroups(7)

fuel.setSigmaT(numpy.array([1.779490E-01, 3.298050E-01, 4.803880E-01, 
                5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01]))
fuel.setSigmaA(numpy.array([8.024800E-03, 3.717400E-03, 2.676900E-02, 
                9.623600E-02, 3.002000E-02, 1.112600E-01, 2.827800E-01]))
fuel.setSigmaS(numpy.array([1.275370E-01, 4.237800E-02, 9.437400E-06, 
                            5.516300E-09, 0., 0., 0., 0., 3.244560E-01,
                            1.631400E-03, 3.142700E-09, 0., 0., 0., 0.,
                            0., 4.509400E-01, 2.679200E-03, 0., 0., 0., 
                            0.,	0., 0., 4.525650E-01, 5.566400E-03, 0., 
                            0., 0., 0., 0., 1.252500E-04, 2.714010E-01,
                            1.025500E-02, 1.002100E-08, 0., 0.,	0., 0.,
                            1.296800E-03, 2.658020E-01, 1.680900E-02, 
                            0., 0., 0., 0., 0.,	8.545800E-03, 2.730800E-01]))

fuel.setSigmaF(numpy.array([7.212060E-03, 8.193010E-04, 6.453200E-03,
                1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01]))
fuel.setNuSigmaF(numpy.array([2.005998E-02, 2.027303E-03, 1.570599E-02, 
                4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01]))
fuel.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 
                1.17610E-07, 0.00000E+00, 0.00000E+00, 0.00000E+00]))

moderator.setSigmaA(numpy.array([6.010500E-04, 1.579300E-05, 3.371600E-04,
                1.940600E-03, 5.741600E-03, 1.500100E-02, 3.723900E-02]))
moderator.setSigmaT(numpy.array([1.592060E-01, 4.129700E-01, 5.903100E-01,
                5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00]))
moderator.setSigmaS(numpy.array([4.447770E-02, 1.134000E-01, 7.234700E-04,
                3.749900E-06, 5.318400E-08, 0., 0., 0., 2.823340E-01,
                1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06,
                1.045500E-06, 0., 0., 3.452560E-01, 2.245700E-01, 
                1.699900E-02, 2.644300E-03, 5.034400E-04, 0., 0., 0.,
                9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02,
                0., 0., 0., 7.143700E-05, 1.391380E-01, 5.118200E-01,
                6.122900E-02, 0., 0., 0., 0., 2.215700E-03, 6.999130E-01,
                5.373200E-01, 0., 0., 0., 0., 0., 1.324400E-01, 2.480700E+00]))
moderator.setSigmaF(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
moderator.setNuSigmaF(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
moderator.setChi(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

log.py_printf('NORMAL', 'Creating surfaces...')

circle1 = Circle(id=1, x=0.0, y=0.0, radius=1.0)
circle2 = Circle(id=2, x=0.0, y=0.0, radius=0.5)
left = XPlane(id=3, x=-2.0)
right = XPlane(id=4, x=2.0)
top = YPlane(id=5, y=2.0)
bottom = YPlane(id=6, y=-2.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(id=1, universe=1, material=2))
cells.append(CellBasic(id=2, universe=1, material=2))
cells.append(CellBasic(id=3, universe=1, material=1))
cells.append(CellFill(id=4, universe=0, universe_fill=2))

cells[0].addSurface(halfspace=-1, surface=circle1)
cells[0].addSurface(halfspace=+1, surface=circle2)
cells[1].addSurface(halfspace=-1, surface=circle2)
cells[2].addSurface(halfspace=+1, surface=circle1)
cells[3].addSurface(halfspace=+1, surface=left)
cells[3].addSurface(halfspace=-1, surface=right)
cells[3].addSurface(halfspace=+1, surface=bottom)
cells[3].addSurface(halfspace=-1, surface=top)

cells[0].setNumRings(3)
cells[1].setNumRings(2)

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=2, width_x=4.0, width_y=4.0)
lattice.setLatticeCells([[1]])

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.addMaterial(fuel)
geometry.addMaterial(moderator)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator()
track_generator.setNumAzim(64)
track_generator.setTrackSpacing(0.05)
track_generator.setGeometry(geometry)
track_generator.generateTracks()

#plotter.plotTracks(track_generator)
plotter.plotSegments(track_generator)
plotter.plotMaterials(geometry, gridsize=250)
plotter.plotCells(geometry, gridsize=250)
plotter.plotFlatSourceRegions(geometry, gridsize=250)

solver = Solver(geometry, track_generator)
solver.setNumThreads(1)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(10000)

