from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from surfaces import *

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from py...')

materials = materialize.materialize('Takeda-materials.py')


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

core_cell = Cell()
core_cell.setFill(materials['core'])

reflector_cell = Cell()
reflector_cell.setFill(materials['reflector'])

control_rod_cell = Cell()
control_rod_cell.setFill(materials['control_rod'])

void_cell = Cell()
void_cell.setFill(materials['void'])

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

core = Universe(name='core cell')
reflector = Universe(name='reflector cell')
control_rod = Universe(name='control rod cell')
void = Universe(name='void cell')
root_universe = Universe(name='root universe')

core.addCell(core_cell)
control_rod.addCell(control_rod_cell)
reflector.addCell(reflector_cell)
void.addCell(void_cell)
root_universe.addCell(root_cell)
