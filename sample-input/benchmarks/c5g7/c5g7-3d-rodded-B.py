import openmoc

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')

openmoc.log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 3D Rodded B' \
                      ' Benchmark Problem...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../../')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

xmin = openmoc.XPlane(x=-32.13, name='xmin')
xmax = openmoc.XPlane(x= 32.13, name='xmax')
ymin = openmoc.YPlane(y=-32.13, name='ymin')
ymax = openmoc.YPlane(y= 32.13, name='ymax')
zmin = openmoc.ZPlane(z=-32.13, name='zmin')
zmax = openmoc.ZPlane(z= 32.13, name='zmax')

xmin.setBoundaryType(openmoc.REFLECTIVE)
xmax.setBoundaryType(openmoc.VACUUM)
ymin.setBoundaryType(openmoc.VACUUM)
ymax.setBoundaryType(openmoc.REFLECTIVE)

# Create ZCylinders for the fuel as well as to discretize the moderator into
# rings
fuel_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.58)
moderator_outer_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.62)


###############################################################################
######################   Creating Cells and Universes   #######################
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

# Control rod pin cell
control_rod_cell = openmoc.Cell()
control_rod_cell.setNumRings(3)
control_rod_cell.setNumSectors(8)
control_rod_cell.setFill(materials['Control Rod'])
control_rod_cell.addSurface(-1, fuel_radius)

control_rod = openmoc.Universe(name='Control Rod')
control_rod.addCell(control_rod_cell)
control_rod.addCell(moderator_ring1)
control_rod.addCell(moderator_ring2)
control_rod.addCell(moderator_ring3)

# Reflector Cells
reflector_cell = openmoc.Cell(name='moderator')
reflector_cell.setFill(materials['Water'])

reflector = openmoc.Universe(name='Reflector')
reflector.addCell(reflector_cell)

refined_reflector_cell = openmoc.Cell(name='Semi-Finely Spaced Reflector')
refined_reflector = openmoc.Universe(name='Semi-Finely Spaced Moderator')
refined_reflector.addCell(refined_reflector_cell)

# Lattice-filled Cells for the assemblies
assembly_uo2_unrod_cell = openmoc.Cell(name='UO2 Assembly Unrodded')
assembly_mox_unrod_cell = openmoc.Cell(name='MOX Assembly Unrodded')
assembly_rfl_unrod_cell_rgt = openmoc.Cell(name='Reflector Unrodded Right')
assembly_rfl_unrod_cell_btm = openmoc.Cell(name='Reflector Unrodded Bottom')
assembly_rfl_unrod_cell_cnr = openmoc.Cell(name='Reflector Unrodded Corner')
assembly_uo2_rod_cell = openmoc.Cell(name='UO2 Assembly Unrodded')
assembly_mox_rod_cell = openmoc.Cell(name='MOX Assembly Unrodded')
assembly_rfl_rod_cell = openmoc.Cell(name='Reflector Rodded')
assembly_rfl_unrod_cell = openmoc.Cell(name='Reflector Unrodded')

# Assembly Universes
assembly_uo2_unrod = openmoc.Universe(name='UO2 Assembly Unrodded')
assembly_mox_unrod = openmoc.Universe(name='MOX Assembly Unrodded')
assembly_rfl_unrod_rgt = openmoc.Universe(name='Rfl Assembly Unrodded Right')
assembly_rfl_unrod_btm = openmoc.Universe(name='Rfl Assembly Unrodded Bottom')
assembly_rfl_unrod_cnr = openmoc.Universe(name='Rfl Assembly Unrodded Corner')
assembly_uo2_rod = openmoc.Universe(name='UO2 Assembly Unrodded')
assembly_mox_rod = openmoc.Universe(name='MOX Assembly Unrodded')
assembly_rfl_rod = openmoc.Universe(name='Rfl Assembly Rodded')
assembly_rfl_unrod = openmoc.Universe(name='Rfl Assembly Unrodded')

assembly_uo2_unrod.addCell(assembly_uo2_unrod_cell)
assembly_mox_unrod.addCell(assembly_mox_unrod_cell)
assembly_rfl_unrod_rgt.addCell(assembly_rfl_unrod_cell_rgt)
assembly_rfl_unrod_btm.addCell(assembly_rfl_unrod_cell_btm)
assembly_rfl_unrod_cnr.addCell(assembly_rfl_unrod_cell_cnr)
assembly_uo2_rod.addCell(assembly_uo2_rod_cell)
assembly_mox_rod.addCell(assembly_mox_rod_cell)
assembly_rfl_rod.addCell(assembly_rfl_rod_cell)
assembly_rfl_unrod.addCell(assembly_rfl_unrod_cell)

# Root Cell/Universe
root_cell = openmoc.Cell(name='Full Geometry')
root_cell.addSurface(halfspace=+1, surface=xmin)
root_cell.addSurface(halfspace=-1, surface=xmax)
root_cell.addSurface(halfspace=+1, surface=ymin)
root_cell.addSurface(halfspace=-1, surface=ymax)
root_cell.addSurface(halfspace=+1, surface=zmin)
root_cell.addSurface(halfspace=-1, surface=zmax)

root_universe = openmoc.Universe(name='Root Universe')
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating lattices...')

lattices = list()

# Sliced up water cells - semi finely spaced
lattices.append(openmoc.Lattice(name='Semi-Finely Spaced Reflector'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
template = [[reflector] * 10] * 10
lattices[-1].setUniverses([template])
refined_reflector_cell.setFill(lattices[-1])


# UO2 unrodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly UO2 Unrodded'))
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
assembly_uo2_unrod_cell.setFill(lattices[-1])


# UO2 rodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly UO2 Rodded'))
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


universes = {1 : uo2, 2 : control_rod, 3 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]

lattices[-1].setUniverses([template])
assembly_uo2_rod_cell.setFill(lattices[-1])


# MOX unrodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly MOX Unrodded'))
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
assembly_mox_unrod_cell.setFill(lattices[-1])


# MOX rodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly MOX Rodded'))
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
             4 : control_rod, 5 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]

lattices[-1].setUniverses([template])
assembly_mox_rod_cell.setFill(lattices[-1])


# Reflector rodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly Reflector Rodded'))
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


universes = {1 : reflector, 2 : control_rod, 3 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]

lattices[-1].setUniverses([template])
assembly_rfl_rod_cell.setFill(lattices[-1])

# Reflector unrodded 17 x 17 assemblies
lattices.append(openmoc.Lattice(name='Assembly Reflector Unrodded'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 17
lattices[-1].setUniverses([template])
assembly_rfl_unrod_cell.setFill(lattices[-1])

# Reflector unrodded assembly right
lattices.append(openmoc.Lattice(name='Assembly Reflector Unrodded Right'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses([template])
assembly_rfl_unrod_cell_rgt.setFill(lattices[-1])

# Reflector unrodded assembly bottom
lattices.append(openmoc.Lattice(name='Assembly Reflector Unrodded Bottom'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
assembly_rfl_unrod_cell_btm.setFill(lattices[-1])

# Reflector unrodded assembly bottom
lattices.append(openmoc.Lattice(name='Assembly Reflector Unrodded Corner'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
assembly_rfl_unrod_cell_cnr.setFill(lattices[-1])


###############################################################################
#########################   Creating Core Lattice   ###########################
###############################################################################

# 3 x 3 x 9 core to represent 3D core
lattices.append(openmoc.Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42, width_z=7.14)
lattices[-1].setUniverses(
  [[[assembly_rfl_rod      , assembly_rfl_rod      , assembly_rfl_unrod_rgt],
    [assembly_rfl_rod      , assembly_rfl_rod      , assembly_rfl_unrod_rgt],
    [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]]]*3+
  [[[assembly_uo2_rod      , assembly_mox_rod      , assembly_rfl_unrod_rgt],
    [assembly_mox_rod      , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
    [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]]]*2+
  [[[assembly_uo2_rod      , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
    [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
    [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]]]*2+
  [[[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
    [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
    [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]]]*2)

# Fill root cell with lattice
root_cell.setFill(lattices[-1])


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.0)
cmfd.setLatticeStructure(51,51)
cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.azim_spacing)
track_generator.setZCoord(-20.0)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry)
openmoc.plotter.plot_cells(geometry)
openmoc.plotter.plot_cmfd_cells(geometry, cmfd)
openmoc.plotter.plot_flat_source_regions(geometry)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
openmoc.plotter.plot_fission_rates(solver)

openmoc.log.py_printf('TITLE', 'Finished')
