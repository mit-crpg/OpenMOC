from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = XPlane(x=-2.0, name='xmin')
xmax = XPlane(x= 2.0, name='xmax')
ymin = YPlane(y=-2.0, name='ymin')
ymax = YPlane(y= 2.0, name='ymax')
zmin = ZPlane(z=-0.5, name='zmin')
zmax = ZPlane(z= 0.5, name='zmax')

xmin.setBoundaryType(PERIODIC)
xmax.setBoundaryType(PERIODIC)
ymin.setBoundaryType(PERIODIC)
ymax.setBoundaryType(PERIODIC)
zmin.setBoundaryType(PERIODIC)
zmax.setBoundaryType(PERIODIC)

large_circle = Circle(x=0.0, y=0.0, radius=0.4, name='large pin')
medium_circle = Circle(x=0.0, y=0.0, radius=0.3, name='medium pin')
small_circle = Circle(x=0.0, y=0.0, radius=0.2, name='small pin')

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

large_fuel = Cell()
large_fuel.setFill(materials['UO2'])
large_fuel.addSurface(halfspace=-1, surface=large_circle)

large_moderator = Cell()
large_moderator.setFill(materials['Water'])
large_moderator.addSurface(halfspace=+1, surface=large_circle)

medium_fuel = Cell()
medium_fuel.setFill(materials['UO2'])
medium_fuel.addSurface(halfspace=-1, surface=medium_circle)

medium_moderator = Cell()
medium_moderator.setFill(materials['Water'])
medium_moderator.addSurface(halfspace=+1, surface=medium_circle)

small_fuel = Cell()
small_fuel.setFill(materials['UO2'])
small_fuel.addSurface(halfspace=-1, surface=small_circle)

small_moderator = Cell()
small_moderator.setFill(materials['Water'])
small_moderator.addSurface(halfspace=+1, surface=small_circle)

root_cell = Cell()
root_cell.addSurface(halfspace=+1, surface=xmin)
root_cell.addSurface(halfspace=-1, surface=xmax)
root_cell.addSurface(halfspace=+1, surface=ymin)
root_cell.addSurface(halfspace=-1, surface=ymax)
root_cell.addSurface(halfspace=+1, surface=zmin)
root_cell.addSurface(halfspace=-1, surface=zmax)

###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

pin1 = Universe(name='large pin cell')
pin2 = Universe(name='medium pin cell')
pin3 = Universe(name='small pin cell')
root_universe = Universe(name='root universe')

pin1.addCell(large_fuel)
pin1.addCell(large_moderator)
pin2.addCell(medium_fuel)
pin2.addCell(medium_moderator)
pin3.addCell(small_fuel)
pin3.addCell(small_moderator)
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 4 x 4 lattice...')

lattice = Lattice(name='4x4 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0, width_z=1.0)
lattice.setUniverses([[pin1, pin2, pin1, pin2],
                      [pin2, pin3, pin2, pin3],
                      [pin1, pin2, pin1, pin2],
                      [pin2, pin3, pin2, pin3]])
root_cell.setFill(lattice)


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setLatticeStructure(2,2,1)
cmfd.setKNearest(1)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()
