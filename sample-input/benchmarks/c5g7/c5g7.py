import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.py_setlevel('NORMAL')

log.py_printf('TITLE', 'Simulating a 4x4 core of 17x17 nested assemblies...')

timer = Timer()

log.py_printf('NORMAL', 'Creating materials...need to use hdf5...')

materials = []
for i in range(7): materials.append(Material(i+1))
for material in materials: material.setNumEnergyGroups(7)

uo2 = materials[0]
uo2.setSigmaT(numpy.array([1.779490E-01, 3.298050E-01, 4.803880E-01, 
                5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01]))
uo2.setSigmaA(numpy.array([8.024800E-03, 3.717400E-03, 2.676900E-02, 
                9.623600E-02, 3.002000E-02, 1.112600E-01, 2.827800E-01]))
uo2.setSigmaS(numpy.array([1.275370E-01, 4.237800E-02, 9.437400E-06, 
                           5.516300E-09, 0., 0., 0., 0., 3.244560E-01,
                           1.631400E-03, 3.142700E-09, 0., 0., 0., 0.,
                           0., 4.509400E-01, 2.679200E-03, 0., 0., 0., 
                           0.,	0., 0., 4.525650E-01, 5.566400E-03, 0., 
                           0., 0., 0., 0., 1.252500E-04, 2.714010E-01,
                           1.025500E-02, 1.002100E-08, 0., 0.,	0., 0.,
                           1.296800E-03, 2.658020E-01, 1.680900E-02, 
                           0., 0., 0., 0., 0., 8.545800E-03, 2.730800E-01]))

uo2.setSigmaF(numpy.array([7.212060E-03, 8.193010E-04, 6.453200E-03,
                1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01]))
uo2.setNuSigmaF(numpy.array([2.005998E-02, 2.027303E-03, 1.570599E-02, 
                4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01]))
uo2.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 
                1.17610E-07, 0., 0., 0.]))

mox43 = materials[1]
mox43.setSigmaA(numpy.array([8.433900E-03, 3.757700E-03, 2.797000E-02, 
                   1.042100E-01, 1.399400E-01, 4.091800E-01, 4.093500E-01]))
mox43.setSigmaT(numpy.array([1.787310E-01, 3.308490E-01, 4.837720E-01, 
                   5.669220E-01, 4.262270E-01, 6.789970E-01, 6.828520E-01]))
mox43.setSigmaS(numpy.array([1.288760E-01, 4.141300E-02, 8.229000E-06, 
                   5.040500E-09, 0., 0., 0., 0., 3.254520E-01, 1.639500E-03, 
                   1.598200E-09, 0., 0., 0., 0., 0., 4.531880E-01, 
                   2.614200E-03, 0., 0., 0., 0., 0., 0., 4.571730E-01, 
                   5.539400E-03, 0., 0., 0., 0., 0., 1.604600E-04, 
                   2.768140E-01, 9.312700E-03, 9.165600E-09, 0., 0., 0.,
                   0., 2.005100E-03, 2.529620E-01, 1.485000E-02, 0., 0.,
                   0., 0., 0., 8.494800E-03, 2.650070E-01]))
mox43.setSigmaF(numpy.array([7.62704E-03, 8.76898E-04, 5.69835E-03, 
                   2.28872E-02, 1.07635E-02, 2.32757E-01, 2.48968E-01]))
mox43.setNuSigmaF(numpy.array([2.175300E-02, 2.535103E-03, 1.626799E-02, 
                   6.547410E-02, 3.072409E-02, 6.666510E-01, 7.139904E-01]))
mox43.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07,
                          0., 0., 0.]))

mox7 = materials[2]
mox7.setSigmaA(numpy.array([9.065700E-03, 4.296700E-03, 3.288100E-02, 
                   1.220300E-01, 1.829800E-01, 5.684600E-01, 5.852100E-01]))
mox7.setSigmaT(numpy.array([1.813230E-01, 3.343680E-01, 4.937850E-01, 
                   5.912160E-01, 4.741980E-01, 8.336010E-01, 8.536030E-01]))
mox7.setSigmaS(numpy.array([1.304570E-01, 4.179200E-02, 8.510500E-06, 
                   5.132900E-09, 0., 0., 0., 0., 3.284280E-01, 1.643600E-03,
                   2.201700E-09, 0., 0., 0., 0., 0., 4.583710E-01, 
                   2.533100E-03, 0., 0., 0., 0., 0., 0., 4.637090E-01,
                   5.476600E-03, 0., 0., 0., 0., 0., 1.761900E-04, 
                   2.823130E-01, 8.728900E-03, 9.001600E-09, 0., 0., 0.,
                   0., 2.276000E-03, 2.497510E-01, 1.311400E-02, 0., 0., 
                   0., 0., 0., 8.864500E-03, 2.595290E-01]))
mox7.setSigmaF(numpy.array([8.25446E-03, 1.32565E-03, 8.42156E-03, 
                   3.28730E-02, 1.59636E-02, 3.23794E-01, 3.62803E-01]))
mox7.setNuSigmaF(numpy.array([2.381395E-02, 3.858689E-03, 2.413400E-02, 
                   9.436622E-02, 4.576988E-02, 9.281814E-01, 1.043200E+00]))
mox7.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 
                    0., 0., 0.]))

mox87 = materials[3]
mox87.setSigmaA(numpy.array([9.486200E-03, 4.655600E-03, 3.624000E-02, 
                    1.327200E-01, 2.084000E-01, 6.587000E-01, 6.901700E-01]))
mox87.setSigmaT(numpy.array([1.830450E-01, 3.367050E-01, 5.005070E-01, 
                    6.061740E-01, 5.027540E-01, 9.210280E-01, 9.552310E-01]))
mox87.setSigmaS(numpy.array([1.315040E-01, 4.204600E-02, 8.697200E-06, 
                    5.193800E-09, 0., 0., 0., 0., 3.304030E-01, 1.646300E-03,
                    2.600600E-09, 0., 0., 0., 0., 0., 4.617920E-01, 
                    2.474900E-03, 0., 0., 0., 0., 0., 0., 4.680210E-01,
                    5.433000E-03, 0., 0., 0., 0., 0., 1.859700E-04,
                    2.857710E-01, 8.397300E-03, 8.928000E-09, 0., 0.,
                    0., 0., 2.391600E-03, 2.476140E-01, 1.232200E-02,
                    0., 0., 0., 0., 0., 8.968100E-03, 2.560930E-01]))
mox87.setSigmaF(numpy.array([8.67209E-03, 1.62426E-03, 1.02716E-02, 
                    3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01]))
mox87.setNuSigmaF(numpy.array([2.518600E-02, 4.739509E-03, 2.947805E-02, 
                    1.122500E-01, 5.530301E-02, 1.074999E+00, 1.239298E+00]))
mox87.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07,
                    0., 0., 0.]))

fiss_chamber = materials[4]
fiss_chamber.setSigmaA(numpy.array([5.113200E-04, 7.580100E-05, 3.157200E-04,
                    1.158200E-03, 3.397500E-03, 9.187800E-03, 2.324200E-02]))
fiss_chamber.setSigmaT(numpy.array([1.260320E-01, 2.931600E-01, 2.842400E-01, 
                    2.809600E-01, 3.344400E-01, 5.656400E-01, 1.172150E+00]))
fiss_chamber.setSigmaS(numpy.array([6.616590E-02, 5.907000E-02, 2.833400E-04, 
                    1.462200E-06, 2.064200E-08, 0., 0., 0., 2.403770E-01,
                    5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06,
                    4.214000E-07,  0., 0., 1.834250E-01, 9.228800E-02, 
                    6.936500E-03, 1.079000E-03, 2.054300E-04, 0., 0., 0., 
                    7.907690E-02, 1.699900E-01, 2.586000E-02, 4.925600E-03,
                    0.,	0., 0., 3.734000E-05, 9.975700E-02, 2.067900E-01, 
                    2.447800E-02, 0., 0., 0., 0., 9.174200E-04, 3.167740E-01, 
                    2.387600E-01, 0., 0., 0., 0., 0., 4.979300E-02, 
                    1.09910E+00]))
fiss_chamber.setSigmaF(numpy.array([4.79002E-09, 5.82564E-09, 4.63719E-07, 
                    5.24406E-06, 1.45390E-07, 7.14972E-07, 2.08041E-06]))
fiss_chamber.setNuSigmaF(numpy.array([1.323401E-08, 1.434500E-08, 1.128599E-06,
                    1.276299E-05, 3.538502E-07, 1.740099E-06, 5.063302E-06]))
fiss_chamber.setChi(numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 
                    1.17610E-07, 0., 0., 0.]))

guide_tube = materials[5]
guide_tube.setSigmaA(numpy.array([5.113200E-04, 7.581300E-05, 3.164300E-04, 
                    1.167500E-03, 3.397700E-03, 9.188600E-03, 2.324400E-02]))
guide_tube.setSigmaT(numpy.array([1.260320E-01, 2.931600E-01, 2.842500E-01, 
                    2.810200E-01, 3.344600E-01, 5.656400E-01, 1.172140E+00]))
guide_tube.setSigmaS(numpy.array([6.616590E-02, 5.907000E-02, 2.833400E-04, 
                    1.462200E-06, 2.064200E-08, 0., 0., 0., 2.403770E-01,
                    5.243500E-02, 2.499000E-04, 1.923900E-05, 2.987500E-06, 
                    4.214000E-07, 0., 0., 1.834250E-01, 9.228800E-02, 
                    6.936500E-03, 1.079000E-03, 2.054300E-04, 0., 0., 0.,
                    7.907690E-02, 1.699900E-01, 2.586000E-02, 4.925600E-03,
                    0., 0., 0., 3.734000E-05, 9.975700E-02, 2.067900E-01,
                    2.447800E-02, 0., 0., 0., 0., 9.174200E-04, 3.167740E-01, 
                    2.387600E-01, 0., 0., 0., 0., 0., 4.979300E-02, 
                    1.099100E+00]))
guide_tube.setSigmaF(numpy.array([0., 0., 0., 0., 0., 0., 0.]))
guide_tube.setNuSigmaF(numpy.array([0., 0., 0., 0., 0., 0., 0.]))
guide_tube.setChi(numpy.array([0., 0., 0., 0., 0., 0., 0.]))

water = materials[6]
water.setSigmaA(numpy.array([6.010500E-04, 1.579300E-05, 3.371600E-04,
                1.940600E-03, 5.741600E-03, 1.500100E-02, 3.723900E-02]))
water.setSigmaT(numpy.array([1.592060E-01, 4.129700E-01, 5.903100E-01,
                5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00]))
water.setSigmaS(numpy.array([4.447770E-02, 1.134000E-01, 7.234700E-04,
                3.749900E-06, 5.318400E-08, 0., 0., 0., 2.823340E-01,
                1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06,
                1.045500E-06, 0., 0., 3.452560E-01, 2.245700E-01, 
                1.699900E-02, 2.644300E-03, 5.034400E-04, 0., 0., 0.,
                9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02,
                0., 0., 0., 7.143700E-05, 1.391380E-01, 5.118200E-01,
                6.122900E-02, 0., 0., 0., 0., 2.215700E-03, 6.999130E-01,
                5.373200E-01, 0., 0., 0., 0., 0., 1.324400E-01, 2.480700E+00]))
water.setSigmaF(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
water.setNuSigmaF(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
water.setChi(numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

                    

log.py_printf('NORMAL', 'Creating surfaces...')

circles = []
planes = []
planes.append(XPlane(id=1, x=-32.13))
planes.append(XPlane(id=2, x=32.13))
planes.append(YPlane(id=3, y=-32.13))
planes.append(YPlane(id=4, y=32.13))
circles.append(Circle(id=5, x=0., y=0., radius=0.54))
circles.append(Circle(id=5, x=0., y=0., radius=0.58))
circles.append(Circle(id=5, x=0., y=0., radius=0.62))
planes[0].setBoundaryType(REFLECTIVE)
planes[1].setBoundaryType(VACUUM)
planes[2].setBoundaryType(VACUUM)
planes[3].setBoundaryType(REFLECTIVE)

log.py_printf('NORMAL', 'Creating cells...')

cells = []
# UO2 pin cells
cells.append(CellBasic(id=1, universe=1, material=1))
cells.append(CellBasic(id=2, universe=1, material=7))
cells.append(CellBasic(id=3, universe=1, material=7))
cells.append(CellBasic(id=4, universe=1, material=7))
cells[0].addSurface(-1, circles[0])
cells[1].addSurface(+1, circles[0])
cells[1].addSurface(-1, circles[1])
cells[2].addSurface(+1, circles[1])
cells[2].addSurface(-1, circles[2])
cells[3].addSurface(+1, circles[2])

# 4.3% MOX pin cells
cells.append(CellBasic(id=5, universe=2, material=2))
cells.append(CellBasic(id=6, universe=2, material=7))
cells.append(CellBasic(id=7, universe=2, material=7))
cells.append(CellBasic(id=8, universe=2, material=7))
cells[4].addSurface(-1, circles[0])
cells[5].addSurface(+1, circles[0])
cells[5].addSurface(-1, circles[1])
cells[6].addSurface(+1, circles[1])
cells[6].addSurface(-1, circles[2])
cells[7].addSurface(+1, circles[2])

# 7% MOX pin cells
cells.append(CellBasic(id=9, universe=3, material=3))
cells.append(CellBasic(id=10, universe=3, material=7))
cells.append(CellBasic(id=11, universe=3, material=7))
cells.append(CellBasic(id=12, universe=3, material=7))
cells[8].addSurface(-1, circles[0])
cells[9].addSurface(+1, circles[0])
cells[9].addSurface(-1, circles[1])
cells[10].addSurface(+1, circles[1])
cells[10].addSurface(-1, circles[2])
cells[11].addSurface(+1, circles[2])

# 8.7% MOX pin cells
cells.append(CellBasic(id=13, universe=4, material=4))
cells.append(CellBasic(id=14, universe=4, material=7))
cells.append(CellBasic(id=15, universe=4, material=7))
cells.append(CellBasic(id=16, universe=4, material=7))
cells[12].addSurface(-1, circles[0])
cells[13].addSurface(+1, circles[0])
cells[13].addSurface(-1, circles[1])
cells[14].addSurface(+1, circles[1])
cells[14].addSurface(-1, circles[2])
cells[15].addSurface(+1, circles[2])

# Fission chamber pin cells
cells.append(CellBasic(id=17, universe=5, material=5))
cells.append(CellBasic(id=18, universe=5, material=7))
cells.append(CellBasic(id=19, universe=5, material=7))
cells.append(CellBasic(id=20, universe=5, material=7))
cells[16].addSurface(-1, circles[0])
cells[17].addSurface(+1, circles[0])
cells[17].addSurface(-1, circles[1])
cells[18].addSurface(+1, circles[1])
cells[18].addSurface(-1, circles[2])
cells[19].addSurface(+1, circles[2])

# Guide tube pin cells
cells.append(CellBasic(id=21, universe=6, material=6))
cells.append(CellBasic(id=22, universe=6, material=7))
cells.append(CellBasic(id=23, universe=6, material=7))
cells.append(CellBasic(id=24, universe=6, material=7))
cells[20].addSurface(-1, circles[0])
cells[21].addSurface(+1, circles[0])
cells[21].addSurface(-1, circles[1])
cells[22].addSurface(+1, circles[1])
cells[22].addSurface(-1, circles[2])
cells[23].addSurface(+1, circles[2])

# Moderator lattice - very finely spaced
cells.append(CellFill(id=25, universe=11, universe_fill=23))

# Moderator lattice - semi-finely spaced
cells.append(CellFill(id=26, universe=12, universe_fill=24))

# Moderator lattice - bottom of geometry
cells.append(CellFill(id=27, universe=13, universe_fill=25))

# Moderator lattice - bottom corner of geometry
cells.append(CellFill(id=28, universe=14, universe_fill=26))

# Moderator lattice right side of geometry
cells.append(CellFill(id=29, universe=15, universe_fill=27))

# Moderator cell
cells.append(CellBasic(id=30, universe=7, material=7))

# Top left, bottom right lattice
cells.append(CellFill(id=32, universe=9, universe_fill=20))

# Top right, bottom left lattice
cells.append(CellFill(id=33, universe=10, universe_fill=21))

# Full geometry
cells.append(CellFill(id=34, universe=0, universe_fill=22))
cells[32].addSurface(+1, planes[0])
cells[32].addSurface(-1, planes[1])
cells[32].addSurface(+1, planes[2])
cells[32].addSurface(-1, planes[3])

log.py_printf('NORMAL', 'Creating assemblies...')

lattices = []

# Top left, bottom right 17 x 17 assemblies
lattices.append(Lattice(id=20, width_x=1.26, width_y=1.26))
lattices[0].setLatticeCells([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
                          [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
                          [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])


# Top right, bottom left 17 x 17 assemblies 
lattices.append(Lattice(id=21, width_x=1.26, width_y=1.26))
lattices[1].setLatticeCells([[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                          [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
                          [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
                          [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
                          [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
                          [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
                          [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                          [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                          [2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2],
                          [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                          [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                          [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
                          [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
                          [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
                          [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
                          [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
                          [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]])

# Sliced up water cells - very finely spaced
lattices.append(Lattice(id=23, width_x=0.063, width_y=0.063))
lattices[2].setLatticeCells([[7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])

# Sliced up water cells - semi finely spaced
lattices.append(Lattice(id=24, width_x=0.126, width_y=0.126))
lattices[3].setLatticeCells([[7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(id=27, width_x=1.26, width_y=1.26))
lattices[4].setLatticeCells([[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7]])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(id=26, width_x=1.26, width_y=1.26))
lattices[5].setLatticeCells([[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(id=25, width_x=1.26, width_y=1.26))
lattices[6].setLatticeCells([[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
                          [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]])


# 4 x 4 core to represent two bundles, water, and vacuum boundary (infinite absorber)
lattices.append(Lattice(id=22, width_x=21.42, width_y=21.42))
lattices[7].setLatticeCells([[9, 10, 15],
                             [10, 9, 15],
                             [13, 13, 14]])


log.py_printf('NORMAL', 'Creating geometry...')

timer.startTimer()
geometry = Geometry()
for material in materials: geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
for lattice in lattices: geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()
timer.stopTimer()
timer.recordSplit('Geometry initialization')
timer.resetTimer()

log.py_printf('NORMAL', 'Initializing the track generator...')

timer.startTimer()
track_generator = TrackGenerator()
track_generator.setNumAzim(256)
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
solver.setNumThreads(16)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(1000)
timer.stopTimer()
timer.recordSplit('Fixed source iteration on host')
timer.printSplits()
