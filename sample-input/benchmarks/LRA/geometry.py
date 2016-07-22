import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('LRA-mgxs.h5', '')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

left = openmoc.XPlane(x=-82.5)
right = openmoc.XPlane(x=82.5)
bottom = openmoc.YPlane(y=-82.5)
top = openmoc.YPlane(y=82.5)
left.setBoundaryType(openmoc.REFLECTIVE)
right.setBoundaryType(openmoc.VACUUM)
bottom.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.VACUUM)


###############################################################################
#                       Creating Cells and Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

# Region 1
region1_cell = openmoc.Cell(name='region 1')
region1_cell.setFill(materials['region_1'])
region1 = openmoc.Universe(name='region 1')
region1.addCell(region1_cell)

# Region 2
region2_cell = openmoc.Cell(name='region 2')
region2_cell.setFill(materials['region_2'])
region2 = openmoc.Universe(name='region 2')
region2.addCell(region2_cell)

# Region 3
region3_cell = openmoc.Cell(name='region 3')
region3_cell.setFill(materials['region_3'])
region3 = openmoc.Universe(name='region 3')
region3.addCell(region3_cell)

# Region 4
region4_cell = openmoc.Cell(name='region 4')
region4_cell.setFill(materials['region_4'])
region4 = openmoc.Universe(name='region 4')
region4.addCell(region4_cell)

# Region 5
region5_cell = openmoc.Cell(name='region 5')
region5_cell.setFill(materials['region_5'])
region5 = openmoc.Universe(name='region 5')
region5.addCell(region5_cell)

# Region 5
region6_cell = openmoc.Cell(name='region 6')
region6_cell.setFill(materials['region_6'])
region6 = openmoc.Universe(name='region 6')
region6.addCell(region6_cell)

# Cells
assembly1_cell = openmoc.Cell(name='assembly 1')
assembly2_cell = openmoc.Cell(name='assembly 2')
assembly3_cell = openmoc.Cell(name='assembly 3')
assembly4_cell = openmoc.Cell(name='assembly 4')
assembly5_cell = openmoc.Cell(name='assembly 5')
assembly6_cell = openmoc.Cell(name='assembly 6')

assembly1 = openmoc.Universe(name='assembly 1')
assembly2 = openmoc.Universe(name='assembly 2')
assembly3 = openmoc.Universe(name='assembly 3')
assembly4 = openmoc.Universe(name='assembly 4')
assembly5 = openmoc.Universe(name='assembly 5')
assembly6 = openmoc.Universe(name='assembly 6')

assembly1.addCell(assembly1_cell)
assembly2.addCell(assembly2_cell)
assembly3.addCell(assembly3_cell)
assembly4.addCell(assembly4_cell)
assembly5.addCell(assembly5_cell)
assembly6.addCell(assembly6_cell)

# Root cell/universe
root_cell = openmoc.Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating LRA lattices...')

# Assembly 1
assembly1_lattice = openmoc.Lattice(name='assembly 1')
assembly1_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region1] * 10] * 10
assembly1_lattice.setUniverses([template])
assembly1_cell.setFill(assembly1_lattice)

# Assembly 2
assembly2_lattice = openmoc.Lattice(name='assembly 2')
assembly2_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region2] * 10] * 10
assembly2_lattice.setUniverses([template])
assembly2_cell.setFill(assembly2_lattice)

# Assembly 3
assembly3_lattice = openmoc.Lattice(name='assembly 3')
assembly3_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region3] * 10] * 10
assembly3_lattice.setUniverses([template])
assembly3_cell.setFill(assembly3_lattice)

# Assembly 4
assembly4_lattice = openmoc.Lattice(name='assembly 4')
assembly4_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region4] * 10] * 10
assembly4_lattice.setUniverses([template])
assembly4_cell.setFill(assembly4_lattice)

# Assembly 5
assembly5_lattice = openmoc.Lattice(name='assembly 5')
assembly5_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region5] * 10] * 10
assembly5_lattice.setUniverses([template])
assembly5_cell.setFill(assembly5_lattice)

# Assembly 6
assembly6_lattice = openmoc.Lattice(name='assembly 6')
assembly6_lattice.setWidth(width_x=1.5, width_y=1.5)
template = [[region6] * 10] * 10
assembly6_lattice.setUniverses([template])
assembly6_cell.setFill(assembly6_lattice)

# Full core
core_lattice = openmoc.Lattice(name='core')
core_lattice.setWidth(width_x=15.0, width_y=15.0)

universes = {7 : assembly1, 8 : assembly2, 9: assembly3,
             10 : assembly4, 11 : assembly5, 12 : assembly6}
template = [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
            [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
            [ 9,  9,  9,  9,  9,  9,  9, 12, 12, 12, 12],
            [ 9,  9,  9,  9,  9,  9,  9, 10, 12, 12, 12],
            [ 8,  7,  7,  7,  7,  8,  8, 11, 11, 12, 12],
            [ 8,  7,  7,  7,  7,  8,  8, 11, 11, 12, 12],
            [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
            [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
            [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
            [ 7,  7,  7,  7,  7,  7,  7,  9,  9, 12, 12],
            [ 8,  7,  7,  7,  7,  8,  8,  9,  9, 12, 12]]

for i in range(11):
  for j in range(11):
    template[i][j] = universes[template[i][j]]
core_lattice.setUniverses([template])
root_cell.setFill(core_lattice)


###############################################################################
#                         Creating the Geometry
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
