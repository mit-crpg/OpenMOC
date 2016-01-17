import openmoc

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

length = 10.0
num_cells_x = 100
num_cells_y = 1

###############################################################################
#######################   Define Material Properties   ########################
###############################################################################

basic_material = openmoc.Material(name='1-group infinite medium')
basic_material.setNumEnergyGroups(1)
basic_material.setSigmaF([0.0414198575])
basic_material.setNuSigmaF([0.0994076580])
basic_material.setSigmaS([0.383259177])
basic_material.setChi([1.0])
basic_material.setSigmaT([0.452648699])

#############################################################################
##########################   Creating Surfaces   ############################
#############################################################################

left = openmoc.XPlane(x=-length/2, name='left')
right = openmoc.XPlane(x=length/2, name='right')
top = openmoc.YPlane(y=length/2, name='top')
bottom = openmoc.YPlane(y=-length/2, name='bottom')

left.setBoundaryType(openmoc.REFLECTIVE)
right.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.REFLECTIVE)

#############################################################################
############################   Creating Cells   #############################
#############################################################################

fill = openmoc.Cell(name='fill')
fill.setFill(basic_material)

root_cell = openmoc.Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)

#############################################################################
##########################    Creating Universes   ##########################
#############################################################################

fill_universe = openmoc.Universe(name='homogeneous fill cell')
fill_universe.addCell(fill)

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(root_cell)

#############################################################################
###########################    Creating Lattices   ##########################
#############################################################################

lattice = openmoc.Lattice(name='MxN lattice')
lattice.setWidth(width_x=length/num_cells_x, width_y=length/num_cells_y)
lattice.setUniverses([[[fill_universe] * num_cells_x]*num_cells_y])
root_cell.setFill(lattice)

#############################################################################
#########################   Creating the Geometry   #########################
#############################################################################

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
