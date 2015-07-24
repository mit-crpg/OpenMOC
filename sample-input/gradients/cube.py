from openmoc import *

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

length = 100.0
num_cells_x = 100 
num_cells_y = 1

###############################################################################
#######################   Define Material Properties   ########################
###############################################################################

basic_material = Material(name='1-group infinite medium')
basic_material.setNumEnergyGroups(1)
basic_material.setSigmaA([0.069389522])
basic_material.setSigmaF([0.0414198575])
basic_material.setNuSigmaF([0.0994076580])
basic_material.setSigmaS([0.383259177])
basic_material.setChi([1.0])
basic_material.setSigmaT([0.452648699])

#############################################################################
##########################   Creating Surfaces   ############################
#############################################################################

left = XPlane(x=-length/2, name='left')
right = XPlane(x=length/2, name='right')
top = YPlane(y=length/2, name='top')
bottom = YPlane(y=-length/2, name='bottom')

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)

#############################################################################
############################   Creating Cells   #############################
#############################################################################

fill = Cell(name='fill')
fill.setFill(basic_material)

root_cell = Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)

#############################################################################
##########################    Creating Universes   ##########################
#############################################################################

fill_universe = Universe(name='homogeneous fill cell')
fill_universe.addCell(fill)

root_universe = Universe(name='root universe')
root_universe.addCell(root_cell)

#############################################################################
###########################    Creating Lattices   ##########################
#############################################################################

lattice = Lattice(name='MxN lattice')
lattice.setWidth(width_x=length/num_cells_x, width_y=length/num_cells_y)
lattice.setUniverses([[fill_universe] * num_cells_x]*num_cells_y)
root_cell.setFill(lattice)

#############################################################################
#########################   Creating the Geometry   #########################
#############################################################################

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
