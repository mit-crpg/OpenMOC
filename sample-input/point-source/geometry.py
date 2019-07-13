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

xmin = openmoc.XPlane(x=-20.0, name='xmin')
xmax = openmoc.XPlane(x= 20.0, name='xmax')
ymin = openmoc.YPlane(y=-20.0, name='ymin')
ymax = openmoc.YPlane(y= 20.0, name='ymax')
zmin = openmoc.ZPlane(z=-20.0, name='zmin')
zmax = openmoc.ZPlane(z= 20.0, name='zmax')

xmin.setBoundaryType(openmoc.VACUUM)
xmax.setBoundaryType(openmoc.VACUUM)
ymin.setBoundaryType(openmoc.VACUUM)
ymax.setBoundaryType(openmoc.VACUUM)
zmin.setBoundaryType(openmoc.VACUUM)
zmax.setBoundaryType(openmoc.VACUUM)


###############################################################################
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

water_cell = openmoc.Cell(name='water')
water_cell.setFill(materials['Water'])

source_cell = openmoc.Cell(name='source')
source_cell.setFill(materials['Water'])

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

water_univ = openmoc.Universe(name='water')
source_univ = openmoc.Universe(name='source')
root_universe = openmoc.Universe(name='root universe')

water_univ.addCell(water_cell)
source_univ.addCell(source_cell)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

# Number of lattice cells
num_x = 100
num_y = 100
num_z = 50

# Compute widths of each lattice cell
width_x = (root_universe.getMaxX() - root_universe.getMinX()) / num_x
width_y = (root_universe.getMaxY() - root_universe.getMinY()) / num_y
width_z = (root_universe.getMaxZ() - root_universe.getMinZ()) / num_z

# Create 2D array of Universes in each lattice cell
universes = [[[water_univ]*num_x for _ in range(num_y)]\
             for _ in range(num_z)]

# Place fixed source Universe at (x=10.0, y=10.0, z=1e-5)
source_x = 10.0
source_y = 10.0
source_z = 1e-5
lat_x = (source_x - root_universe.getMinX()) / width_x
lat_y = (source_y - root_universe.getMinY()) / width_y
lat_z = (root_universe.getMaxZ() - source_z) / width_z
universes[int(lat_z)][int(lat_y)][int(lat_x)] = source_univ

log.py_printf('NORMAL', 'Creating a {0}x{1}x{2} lattice...'.\
              format(num_x, num_y, num_z))

lattice = openmoc.Lattice(name='{0}x{1}x{2} lattice'.format(num_x, num_y, num_z))
lattice.setWidth(width_x=width_x, width_y=width_y, width_z=width_z)
lattice.setUniverses(universes)
root_cell.setFill(lattice)


###############################################################################
#                         Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
