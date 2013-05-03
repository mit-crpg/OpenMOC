import numpy
from openmoc import *
import openmoc.log as log
import openmoc.materialize as materialize
import openmoc.process as process


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

num_threads = 4
track_spacing = 0.1
num_azim = 16
tolerance = 1E-3
max_iters = 1000
gridsize = 500

setOutputDirectory('Tiny-Lattice-Strong-Scaling')

log.setLogLevel('INFO')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.hdf5')

uo2_id = materials['UO2'].getId()
water_id = materials['Water'].getId()


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=0.8)
left = XPlane(x=-2.0)
right = XPlane(x=2.0)
top = YPlane(y=2.0)
bottom = YPlane(y=-2.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=uo2_id))
cells.append(CellBasic(universe=1, material=water_id))
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

log.py_printf('NORMAL', 'Creating simple 2 x 2 lattice...')

lattice = Lattice(id=2, width_x=2.0, width_y=2.0)
lattice.setLatticeCells([[1, 1], [1, 1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)

geometry.initializeFlatSourceRegions()


###############################################################################
##############################   Strong Scaling   #############################
###############################################################################

log.py_printf('NORMAL', 'Running a strong scaling study...')

process.strongScalingStudy(geometry, precision='all', compiler='all',
                           num_threads=numpy.linspace(1,24,24))
process.weakScalingStudy(geometry, precision='all', compiler='all',
                           num_threads=numpy.linspace(1,24,24))
