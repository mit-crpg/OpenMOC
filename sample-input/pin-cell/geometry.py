from mpi4py import MPI
import openmoc
openmoc.log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
xmin = openmoc.XPlane(x=-2.0, name='xmin')
ymin = openmoc.YPlane(y=-2.0, name='ymin')
zmin = openmoc.ZPlane(z=-2.0, name='zmin')
xmax = openmoc.XPlane(x=2.0, name='xmax')
ymax = openmoc.YPlane(y=2.0, name='ymax')
zmax = openmoc.ZPlane(z=2.0, name='zmax')

xmin.setBoundaryType(openmoc.PERIODIC)
ymin.setBoundaryType(openmoc.PERIODIC)
zmin.setBoundaryType(openmoc.PERIODIC)
xmax.setBoundaryType(openmoc.PERIODIC)
ymax.setBoundaryType(openmoc.PERIODIC)
zmax.setBoundaryType(openmoc.PERIODIC)

###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

fuel = openmoc.Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=zcylinder)
fuel.addSurface(halfspace=+1, surface=zmin)
fuel.addSurface(halfspace=-1, surface=zmax)

moderator = openmoc.Cell(name='moderator')
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

openmoc.log.py_printf('NORMAL', 'Creating universes...')

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(fuel)
root_universe.addCell(moderator)


###############################################################################
#                         Creating the Geometry
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
print dir(geometry)
geometry.setRootUniverse(root_universe)
geometry.setDomainDecomposition(1, 1, 2)
geometry.initializeFlatSourceRegions()
