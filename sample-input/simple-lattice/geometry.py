import openmoc
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = openmoc.XPlane(x=-2.0, name='xmin')
xmax = openmoc.XPlane(x= 2.0, name='xmax')
ymin = openmoc.YPlane(y=-2.0, name='ymin')
ymax = openmoc.YPlane(y= 2.0, name='ymax')
zmin = openmoc.ZPlane(z=-10.0, name='zmin')
zmax = openmoc.ZPlane(z= 10.0, name='zmax')

xmin.setBoundaryType(openmoc.REFLECTIVE)
xmax.setBoundaryType(openmoc.REFLECTIVE)
ymin.setBoundaryType(openmoc.REFLECTIVE)
ymax.setBoundaryType(openmoc.REFLECTIVE)
zmin.setBoundaryType(openmoc.REFLECTIVE)
zmax.setBoundaryType(openmoc.REFLECTIVE)

large_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.4, name='large pin')
medium_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.3, name='medium pin')
small_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.2, name='small pin')

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

large_fuel = openmoc.Cell()
large_fuel.setFill(materials['UO2'])
large_fuel.addSurface(halfspace=-1, surface=large_zcylinder)

large_moderator = openmoc.Cell()
large_moderator.setFill(materials['Water'])
large_moderator.addSurface(halfspace=+1, surface=large_zcylinder)

medium_fuel = openmoc.Cell()
medium_fuel.setFill(materials['UO2'])
medium_fuel.addSurface(halfspace=-1, surface=medium_zcylinder)

medium_moderator = openmoc.Cell()
medium_moderator.setFill(materials['Water'])
medium_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)

small_fuel = openmoc.Cell()
small_fuel.setFill(materials['UO2'])
small_fuel.addSurface(halfspace=-1, surface=small_zcylinder)

small_moderator = openmoc.Cell()
small_moderator.setFill(materials['Water'])
small_moderator.addSurface(halfspace=+1, surface=small_zcylinder)

root_cell = openmoc.Cell()
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

pin1 = openmoc.Universe(name='large pin cell')
pin2 = openmoc.Universe(name='medium pin cell')
pin3 = openmoc.Universe(name='small pin cell')
root_universe = openmoc.Universe(name='root universe')

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

lattice = openmoc.Lattice(name='4x4 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0, width_z=20.0)
lattice.setUniverses([[[pin1, pin2, pin1, pin2],
                       [pin2, pin3, pin2, pin3],
                       [pin1, pin2, pin1, pin2],
                       [pin2, pin3, pin2, pin3]]])
root_cell.setFill(lattice)


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setLatticeStructure(2,2,4)
cmfd.setCMFDRelaxationFactor(0.7)
cmfd.setKNearest(3)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()
