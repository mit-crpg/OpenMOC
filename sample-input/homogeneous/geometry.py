import openmoc
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

log.set_log_level('NORMAL')

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = openmoc.XPlane(x=-5.0, name='xmin')
xmax = openmoc.XPlane(x= 5.0, name='xmax')
ymin = openmoc.YPlane(y=-5.0, name='ymin')
ymax = openmoc.YPlane(y= 5.0, name='ymax')
zmin = openmoc.ZPlane(z=-5.0, name='zmin')
zmax = openmoc.ZPlane(z= 5.0, name='zmax')

xmin.setBoundaryType(openmoc.REFLECTIVE)
xmax.setBoundaryType(openmoc.REFLECTIVE)
ymin.setBoundaryType(openmoc.REFLECTIVE)
ymax.setBoundaryType(openmoc.REFLECTIVE)
zmin.setBoundaryType(openmoc.REFLECTIVE)
zmax.setBoundaryType(openmoc.REFLECTIVE)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = openmoc.Cell(name='fuel')
fuel.setFill(materials['UO2'])

moderator = openmoc.Cell(name='moderator')
moderator.setFill(materials['UO2'])

root_cell = openmoc.Cell(name='root cell')
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

fue_univ = openmoc.Universe(name='homogeneous fue cell')
fue_univ.addCell(fuel)

mod_univ = openmoc.Universe(name='homogeneous mod cell')
mod_univ.addCell(moderator)

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 10 x 10 lattice...')

f = fue_univ

lattice = openmoc.Lattice(name='10x10 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0, width_z=1.0)
lattice.setUniverses([[[f, f, f, f, f, f, f, f, f, f],
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

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
