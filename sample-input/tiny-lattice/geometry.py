from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter

log.set_log_level('NORMAL')


###############################################################################
#                              Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                              Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

zcylinder = ZCylinder(x=0.0, y=0.0, radius=0.8, name='pin')
xmin = XPlane(x=-2.0, name='xmin')
xmax = XPlane(x= 2.0, name='xmax')
ymin = YPlane(y=-2.0, name='ymin')
ymax = YPlane(y= 2.0, name='ymax')
zmin = ZPlane(z=-1.0, name='zmin')
zmax = ZPlane(z= 1.0, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)


###############################################################################
#                                Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=zcylinder)

moderator = Cell(name='moderator')
moderator.setFill(materials['Water'])
moderator.addSurface(halfspace=+1, surface=zcylinder)

root_cell = Cell(name='root cell')
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

pincell = Universe(name='pin cell')
root_universe = Universe(name='root universe')

pincell.addCell(fuel)
pincell.addCell(moderator)
root_universe.addCell(root_cell)


###############################################################################
#                             Creating Lattices
###############################################################################

log.py_printf('NORMAL', 'Creating simple 2 x 2 lattice...')

lattice = Lattice(name='2x2 lattice')
lattice.setWidth(width_x=2.0, width_y=2.0, width_z=2.0)
lattice.setUniverses([[[pincell, pincell], [pincell, pincell]]])

root_cell.setFill(lattice)


###############################################################################
#                            Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
