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

xmin = XPlane(x=-5.0, name='xmin')
xmax = XPlane(x= 5.0, name='xmax')
ymin = YPlane(y=-5.0, name='ymin')
ymax = YPlane(y= 5.0, name='ymax')
zmin = ZPlane(z=-5.0, name='zmin')
zmax = ZPlane(z= 5.0, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = Cell(name='fuel')
fuel.setFill(materials['UO2'])

moderator = Cell(name='moderator')
moderator.setFill(materials['UO2'])

root_cell = Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=xmin)
root_cell.addSurface(halfspace=-1, surface=xmax)
root_cell.addSurface(halfspace=+1, surface=ymin)
root_cell.addSurface(halfspace=-1, surface=ymax)
root_cell.addSurface(halfspace=+1, surface=zmin)
root_cell.addSurface(halfspace=-1, surface=zmax)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

fue_univ = Universe(name='homogeneous fue cell')
fue_univ.addCell(fuel)

mod_univ = Universe(name='homogeneous mod cell')
mod_univ.addCell(moderator)

root_universe = Universe(name='root universe')
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 10 x 10 lattice...')

f = fue_univ

lattice = Lattice(name='10x10 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0, width_z=1.0)
lattice.setUniverses3D([[[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]],
                        [[f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f],
                         [f, f, f, f, f, f, f, f, f, f]]])

root_cell.setFill(lattice)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
