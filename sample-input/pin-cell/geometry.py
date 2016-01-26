from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.py')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

zcylinder = ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
xmin = XPlane(x=-2.0, name='xmin')
ymin = YPlane(y=-2.0, name='ymin')
zmin = ZPlane(z=-2.0, name='zmin')
xmax = XPlane(x=2.0, name='xmax')
ymax = YPlane(y=2.0, name='ymax')
zmax = ZPlane(z=2.0, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

###############################################################################
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=zcylinder)
fuel.addSurface(halfspace=+1, surface=zmin)
fuel.addSurface(halfspace=-1, surface=zmax)

moderator = Cell(name='moderator')
moderator.setFill(materials['Water'])
moderator.addSurface(halfspace=+1, surface=zcylinder)
moderator.addSurface(halfspace=+1, surface=xmin)
moderator.addSurface(halfspace=-1, surface=xmax)
moderator.addSurface(halfspace=+1, surface=ymin)
moderator.addSurface(halfspace=-1, surface=ymax)
moderator.addSurface(halfspace=+1, surface=zmin)
moderator.addSurface(halfspace=-1, surface=zmax)

###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(fuel)
root_universe.addCell(moderator)


###############################################################################
#                         Creating the Geometry
###############################################################################


log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
