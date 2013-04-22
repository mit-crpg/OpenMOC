import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('INFO')

log.py_printf('TITLE', 'Simulating a 4x4 core of 17x17 nested assemblies...')

timer = Timer()

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

circles = []
planes = []
radii = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
planes.append(XPlane(id=1, x=-34.0))
planes.append(XPlane(id=2, x=34.0))
planes.append(YPlane(id=3, y=-34.0))
planes.append(YPlane(id=4, y=34.0))
for plane in planes: plane.setBoundaryType(REFLECTIVE)

cell_id = 5
for r in radii: 
    circles.append(Circle(id=cell_id, x=0.0, y=0.0, radius=r))
    cell_id += 1

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(id=1, universe=1, material=1))
cells.append(CellBasic(id=2, universe=1, material=2))
cells.append(CellBasic(id=3, universe=2, material=1))
cells.append(CellBasic(id=4, universe=2, material=2))
cells.append(CellBasic(id=5, universe=3, material=1))
cells.append(CellBasic(id=6, universe=3, material=2))
cells.append(CellBasic(id=7, universe=4, material=1))
cells.append(CellBasic(id=8, universe=4, material=2))
cells.append(CellBasic(id=9, universe=5, material=1))
cells.append(CellBasic(id=10, universe=5, material=2))
cells.append(CellBasic(id=11, universe=6, material=1))
cells.append(CellBasic(id=12, universe=6, material=2))

cells.append(CellFill(id=13, universe=9, universe_fill=7))
cells.append(CellFill(id=14, universe=11, universe_fill=8))
cells.append(CellFill(id=15, universe=0, universe_fill=10))

cells[0].addSurface(halfspace=-1, surface=circles[0])
cells[1].addSurface(halfspace=+1, surface=circles[0])
cells[2].addSurface(halfspace=-1, surface=circles[1])
cells[3].addSurface(halfspace=+1, surface=circles[1])
cells[4].addSurface(halfspace=-1, surface=circles[2])
cells[5].addSurface(halfspace=+1, surface=circles[2])
cells[6].addSurface(halfspace=-1, surface=circles[3])
cells[7].addSurface(halfspace=+1, surface=circles[3])
cells[8].addSurface(halfspace=-1, surface=circles[4])
cells[9].addSurface(halfspace=+1, surface=circles[4])
cells[10].addSurface(halfspace=-1, surface=circles[5])
cells[11].addSurface(halfspace=+1, surface=circles[5])

cells[14].addSurface(halfspace=+1, surface=planes[0])
cells[14].addSurface(halfspace=-1, surface=planes[1])
cells[14].addSurface(halfspace=+1, surface=planes[2])
cells[14].addSurface(halfspace=-1, surface=planes[3])

log.py_printf('NORMAL', 'Creating nested lattices...')

# 1st 17x17 assembly
assembly1 = Lattice(id=7, width_x=1.0, width_y=1.0)
assembly1.setLatticeCells([[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],
                          [2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2],
                          [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1]])

# 2nd 17x17 assembly
assembly2 = Lattice(id=8, width_x=1.0, width_y=1.0)
assembly2.setLatticeCells([[4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
                          [5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5],
                          [4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4]])

# 4x4 core
core = Lattice(id=10, width_x=17.0, width_y=17.0)
core.setLatticeCells([[9, 11, 9, 11],
                      [11, 9, 11, 9],
                      [9, 11, 9, 11],
                      [11, 9, 11, 9]])

log.py_printf('NORMAL', 'Creating geometry...')

timer.startTimer()
geometry = Geometry()
geometry.addMaterial(fuel)
geometry.addMaterial(moderator)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(core)
geometry.addLattice(assembly1)
geometry.addLattice(assembly2)
geometry.initializeFlatSourceRegions()
timer.stopTimer()
timer.recordSplit('Geometry initialization')
timer.resetTimer()

log.py_printf('NORMAL', 'Initializing the track generator...')

timer.startTimer()
track_generator = TrackGenerator()
track_generator.setNumAzim(16)
track_generator.setTrackSpacing(0.1)
track_generator.setGeometry(geometry)
track_generator.generateTracks()
timer.stopTimer()
timer.recordSplit('Generating tracks')
timer.resetTimer()

#plotter.plotTracks(track_generator)
#plotter.plotMaterials(geometry, gridsize=250)
#plotter.plotCells(geometry, gridsize=250)
#plotter.plotFlatSourceRegions(geometry, gridsize=250)

timer.startTimer()
solver = Solver(geometry, track_generator)
solver.setNumThreads(2)
solver.setSourceConvergenceThreshold(1E-4)
solver.convergeSource(1000)
timer.stopTimer()
timer.recordSplit('Fixed source iteration on host')
timer.printSplits()
