from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

log.set_log_level('NORMAL')
refines = 1

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from py...')

materials = materialize.materialize('Takeda-materials.py')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = XPlane(x=-12.5, name='xmin')
xmax = XPlane(x= 12.5, name='xmax')
ymin = YPlane(y=-12.5, name='ymin')
ymax = YPlane(y= 12.5, name='ymax')
zmin = ZPlane(z=-12.5, name='zmin')
zmax = ZPlane(z= 12.5, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(VACUUM)
ymax.setBoundaryType(VACUUM)
zmax.setBoundaryType(VACUUM)


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


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 4 x 4 lattice...')

c = core
v = void
a = control_rod
r = reflector

lattice = Lattice(name='4x4 lattice')
lattice.setWidth(width_x=5.0/refines, width_y=5.0/refines, width_z=5.0/refines)
template = [[[r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, v, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, v, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, v, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, v, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, v, r]]]

# refine lattice
template_refined = [[ []*5 ]*5 ]*5
for k in range(5):
  for j in range(5):
    a = template[k][j]
    template_refined[k][j] = sum([[a[i]]*refines for i in range(len(a))], [])

lattice.setUniverses3D(template_refined)
root_cell.setFill(lattice)


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setMOCRelaxationFactor(1.0)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(5, 5, 5)
cmfd.setKNearest(4)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()
