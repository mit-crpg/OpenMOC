import openmoc

###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../../')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

left = openmoc.XPlane(x=-32.13, name='left')
right = openmoc.XPlane(x=32.13, name='right')
top = openmoc.YPlane(y=32.13, name='top')
bottom = openmoc.YPlane(y=-32.13, name='bottom')
left.setBoundaryType(openmoc.REFLECTIVE)
right.setBoundaryType(openmoc.VACUUM)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.VACUUM)
boundaries = [left, right, top, bottom]

# Create ZCylinders for the fuel as well as to discretize the moderator into
# rings
fuel_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.58)
moderator_outer_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.62)


###############################################################################
#                        Creating Cells and Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator_ring1 = openmoc.Cell()
moderator_ring2 = openmoc.Cell()
moderator_ring3 = openmoc.Cell()
moderator_ring1.setNumSectors(8)
moderator_ring2.setNumSectors(8)
moderator_ring3.setNumSectors(8)
moderator_ring1.setFill(materials['Water'])
moderator_ring2.setFill(materials['Water'])
moderator_ring3.setFill(materials['Water'])
moderator_ring1.addSurface(+1, fuel_radius)
moderator_ring1.addSurface(-1, moderator_inner_radius)
moderator_ring2.addSurface(+1, moderator_inner_radius)
moderator_ring2.addSurface(-1, moderator_outer_radius)
moderator_ring3.addSurface(+1, moderator_outer_radius)

# UO2 pin cell
uo2_cell = openmoc.Cell()
uo2_cell.setNumRings(3)
uo2_cell.setNumSectors(8)
uo2_cell.setFill(materials['UO2'])
uo2_cell.addSurface(-1, fuel_radius)

uo2 = openmoc.Universe(name='UO2')
uo2.addCell(uo2_cell)
uo2.addCell(moderator_ring1)
uo2.addCell(moderator_ring2)
uo2.addCell(moderator_ring3)

# 4.3% MOX pin cell
mox43_cell = openmoc.Cell()
mox43_cell.setNumRings(3)
mox43_cell.setNumSectors(8)
mox43_cell.setFill(materials['MOX-4.3%'])
mox43_cell.addSurface(-1, fuel_radius)

mox43 = openmoc.Universe(name='MOX-4.3%')
mox43.addCell(mox43_cell)
mox43.addCell(moderator_ring1)
mox43.addCell(moderator_ring2)
mox43.addCell(moderator_ring3)

# 7% MOX pin cell
mox7_cell = openmoc.Cell()
mox7_cell.setNumRings(3)
mox7_cell.setNumSectors(8)
mox7_cell.setFill(materials['MOX-7%'])
mox7_cell.addSurface(-1, fuel_radius)

mox7 = openmoc.Universe(name='MOX-7%')
mox7.addCell(mox7_cell)
mox7.addCell(moderator_ring1)
mox7.addCell(moderator_ring2)
mox7.addCell(moderator_ring3)

# 8.7% MOX pin cell
mox87_cell = openmoc.Cell()
mox87_cell.setNumRings(3)
mox87_cell.setNumSectors(8)
mox87_cell.setFill(materials['MOX-8.7%'])
mox87_cell.addSurface(-1, fuel_radius)

mox87 = openmoc.Universe(name='MOX-8.7%')
mox87.addCell(mox87_cell)
mox87.addCell(moderator_ring1)
mox87.addCell(moderator_ring2)
mox87.addCell(moderator_ring3)

# Fission chamber pin cell
fission_chamber_cell = openmoc.Cell()
fission_chamber_cell.setNumRings(3)
fission_chamber_cell.setNumSectors(8)
fission_chamber_cell.setFill(materials['Fission Chamber'])
fission_chamber_cell.addSurface(-1, fuel_radius)

fission_chamber = openmoc.Universe(name='Fission Chamber')
fission_chamber.addCell(fission_chamber_cell)
fission_chamber.addCell(moderator_ring1)
fission_chamber.addCell(moderator_ring2)
fission_chamber.addCell(moderator_ring3)

# Guide tube pin cell
guide_tube_cell = openmoc.Cell()
guide_tube_cell.setNumRings(3)
guide_tube_cell.setNumSectors(8)
guide_tube_cell.setFill(materials['Guide Tube'])
guide_tube_cell.addSurface(-1, fuel_radius)

guide_tube = openmoc.Universe(name='Guide Tube')
guide_tube.addCell(guide_tube_cell)
guide_tube.addCell(moderator_ring1)
guide_tube.addCell(moderator_ring2)
guide_tube.addCell(moderator_ring3)

# Reflector
reflector_cell = openmoc.Cell(name='moderator')
reflector_cell.setFill(materials['Water'])

reflector = openmoc.Universe(name='Reflector')
reflector.addCell(reflector_cell)

# Cells
assembly1_cell = openmoc.Cell(name='Assembly 1')
assembly2_cell = openmoc.Cell(name='Assembly 2')
refined_reflector_cell = openmoc.Cell(name='Semi-Finely Spaced Reflector')
right_reflector_cell = openmoc.Cell(name='Right Reflector')
corner_reflector_cell = openmoc.Cell(name='Bottom Corner Reflector')
bottom_reflector_cell = openmoc.Cell(name='Bottom Reflector')

assembly1 = openmoc.Universe(name='Assembly 1')
assembly2 = openmoc.Universe(name='Assembly 2')
refined_reflector = openmoc.Universe(name='Semi-Finely Spaced Moderator')
right_reflector = openmoc.Universe(name='Right Reflector')
corner_reflector = openmoc.Universe(name='Bottom Corner Reflector')
bottom_reflector = openmoc.Universe(name='Bottom Reflector')

assembly1.addCell(assembly1_cell)
assembly2.addCell(assembly2_cell)
refined_reflector.addCell(refined_reflector_cell)
right_reflector.addCell(right_reflector_cell)
corner_reflector.addCell(corner_reflector_cell)
bottom_reflector.addCell(bottom_reflector_cell)

# Root Cell/Universe
root_cell = openmoc.Cell(name='Full Geometry')
root_cell.addSurface(+1, boundaries[0])
root_cell.addSurface(-1, boundaries[1])
root_cell.addSurface(-1, boundaries[2])
root_cell.addSurface(+1, boundaries[3])

root_universe = openmoc.Universe(name='Root Universe')
root_universe.addCell(root_cell)


###############################################################################
#                             Creating Lattices
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating lattices...')

lattices = list()

# Top left, bottom right 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly 1'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1],
            [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
            [1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

universes = {1 : uo2, 2 : guide_tube, 3 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses([template])
assembly1_cell.setFill(lattices[-1])

# Top right, bottom left 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly 2'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
universes = {1 : mox43, 2 : mox7, 3 : mox87,
             4 : guide_tube, 5 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses([template])
assembly2_cell.setFill(lattices[-1])

# Sliced up water cells - semi finely spaced
lattices.append(openmoc.Lattice(name='Semi-Finely Spaced Reflector'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
template = [[reflector] * 10] * 10
lattices[-1].setUniverses([template])
refined_reflector_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(openmoc.Lattice(name='Right Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses([template])
right_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(openmoc.Lattice(name='Bottom Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom of geometry
lattices.append(openmoc.Lattice(name='Bottom Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
bottom_reflector_cell.setFill(lattices[-1])

# 4 x 4 core to represent two bundles and water
lattices.append(openmoc.Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)
lattices[-1].setUniverses([[
     [assembly1,        assembly2,        right_reflector],
     [assembly2,        assembly1,        right_reflector],
     [bottom_reflector, bottom_reflector, corner_reflector]]])
root_cell.setFill(lattices[-1])


###############################################################################
#                         Creating the Geometry
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()
