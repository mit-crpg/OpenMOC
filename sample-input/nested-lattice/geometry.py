import openmoc
import openmoc.log as log
import openmoc.plotter as plotter

log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = openmoc.XPlane(x=-2.0, name='xmin')
xmax = openmoc.XPlane(x= 2.0, name='xmax')
ymin = openmoc.YPlane(y=-2.0, name='ymin')
ymax = openmoc.YPlane(y= 2.0, name='ymax')
zmin = openmoc.ZPlane(z=-0.5, name='zmin')
zmax = openmoc.ZPlane(z= 0.5, name='zmax')

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
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

large_fuel = openmoc.Cell(name='large pin fuel')
large_fuel.setNumRings(3)
large_fuel.setNumSectors(8)
large_fuel.setFill(materials['UO2'])
large_fuel.addSurface(halfspace=-1, surface=large_zcylinder)

large_moderator = openmoc.Cell(name='large pin moderator')
large_moderator.setNumSectors(8)
large_moderator.setFill(materials['Water'])
large_moderator.addSurface(halfspace=+1, surface=large_zcylinder)

medium_fuel = openmoc.Cell(name='medium pin fuel')
medium_fuel.setNumRings(3)
medium_fuel.setNumSectors(8)
medium_fuel.setFill(materials['UO2'])
medium_fuel.addSurface(halfspace=-1, surface=medium_zcylinder)

medium_moderator = openmoc.Cell(name='medium pin moderator')
medium_moderator.setNumSectors(8)
medium_moderator.setFill(materials['Water'])
medium_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)

small_fuel = openmoc.Cell(name='small pin fuel')
small_fuel.setNumRings(3)
small_fuel.setNumSectors(8)
small_fuel.setFill(materials['UO2'])
small_fuel.addSurface(halfspace=-1, surface=small_zcylinder)

small_moderator = openmoc.Cell(name='small pin moderator')
small_moderator.setNumSectors(8)
small_moderator.setFill(materials['Water'])
small_moderator.addSurface(halfspace=+1, surface=small_zcylinder)

lattice_cell = openmoc.Cell(name='lattice cell')

root_cell = openmoc.Cell(name='root cell')
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
assembly = openmoc.Universe(name='2x2 lattice')
root_universe = openmoc.Universe(name='root universe')

pin1.addCell(large_fuel)
pin1.addCell(large_moderator)
pin2.addCell(medium_fuel)
pin2.addCell(medium_moderator)
pin3.addCell(small_fuel)
pin3.addCell(small_moderator)
assembly.addCell(lattice_cell)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

log.py_printf('NORMAL', 'Creating nested 2 x 2 lattices...')

# 2x2 assembly
lattice = openmoc.Lattice(name='2x2 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0)
lattice.setUniverses([[[pin1, pin2], [pin1, pin3]]])
lattice_cell.setFill(lattice)

# 2x2 core
core = openmoc.Lattice(name='2x2 core')
core.setWidth(width_x=2.0, width_y=2.0, width_z=1.0)
core.setUniverses([[[assembly, assembly], [assembly, assembly]]])
root_cell.setFill(core)


###############################################################################
#                         Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
