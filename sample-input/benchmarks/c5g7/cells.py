import openmoc
import openmoc.materialize as materialize
from surfaces import surfaces

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

materials = materialize.materialize('../../c5g7-materials.py')

###############################################################################
######################   Creating Cells and Universes   #######################
###############################################################################

cells = {}

# Instantiate Cells
cells['UO2']             = opemoc.Cell()
cells['MOX 4.3%']        = opemoc.Cell()
cells['MOX 7.0%']        = opemoc.Cell()
cells['MOX 8.7%']        = opemoc.Cell()
cells['Guide Tube']      = opemoc.Cell()
cells['Fission Chamber'] = opemoc.Cell()
cells['Control Rod']     = opemoc.Cell()
cells['Moderator']       = opemoc.Cell()
cells['Reflector']       = opemoc.Cell()

# Set material fills
cells['UO2']            .setFill(materials['UO2'])
cells['MOX 4.3%']       .setFill(materials['MOX 4.3%'])
cells['MOX 7.0%']       .setFill(materials['MOX 7.0%'])
cells['MOX 8.7%']       .setFill(materials['MOX 8.7%'])
cells['Guide Tube']     .setFill(materials['Guide Tube'])
cells['Fission Chamber'].setFill(materials['Fission Chamber'])
cells['Control Rod']    .setFill(materials['Control Rod'])
cells['Moderator']      .setFill(materials['Water'])
cells['Reflector']      .setFill(materials['Water'])

# Set rings and sectors
cells['UO2']            .setNumRings(3)
cells['MOX 4.3%']       .setNumRings(3)
cells['MOX 7.0%']       .setNumRings(3)
cells['MOX 8.7%']       .setNumRings(3)
cells['Guide Tube']     .setNumRings(3)
cells['Fission Chamber'].setNumRings(3)
cells['Control Rod']    .setNumRings(3)
cells['Moderator']      .setNumRings(2)
cells['UO2']            .setNumSectors(8)
cells['MOX 4.3%']       .setNumSectors(8)
cells['MOX 7.0%']       .setNumSectors(8)
cells['MOX 8.7%']       .setNumSectors(8)
cells['Guide Tube']     .setNumSectors(8)
cells['Fission Chamber'].setNumSectors(8)
cells['Control Rod']    .setNumSectors(8)
cells['Moderator']      .setNumSectors(8)

# Add surfaces
cells['UO2']            .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 4.3%']       .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 7.0%']       .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 8.7%']       .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Guide Tube']     .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Fission Chamber'].addSurface(-1, surfaces['Fuel Cylinder'])
cells['Control Rod']    .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Moderator']      .addSurface(+1, surfaces['Fuel Cylinder'])

reflector = Universe(name='Reflector')
reflector.addCell(reflector_cell)

refined_reflector_cell = Cell(name='Semi-Finely Spaced Reflector')
refined_reflector = Universe(name='Semi-Finely Spaced Moderator')
refined_reflector.addCell(refined_reflector_cell)

# Cells
assembly_uo2_unrod_cell = Cell(name='UO2 Assembly Unrodded')
assembly_mox_unrod_cell = Cell(name='MOX Assembly Unrodded')
assembly_rfl_unrod_cell_rgt = Cell(name='Reflector Unrodded Right')
assembly_rfl_unrod_cell_btm = Cell(name='Reflector Unrodded Bottom')
assembly_rfl_unrod_cell_cnr = Cell(name='Reflector Unrodded Corner')
assembly_uo2_rod_cell = Cell(name='UO2 Assembly Unrodded')
assembly_mox_rod_cell = Cell(name='MOX Assembly Unrodded')
assembly_rfl_rod_cell = Cell(name='Reflector Rodded')
assembly_rfl_unrod_cell = Cell(name='Reflector Unrodded')

assembly_uo2_unrod = Universe(name='UO2 Assembly Unrodded')
assembly_mox_unrod = Universe(name='MOX Assembly Unrodded')
assembly_rfl_unrod_rgt = Universe(name='Rfl Assembly Unrodded Right')
assembly_rfl_unrod_btm = Universe(name='Rfl Assembly Unrodded Bottom')
assembly_rfl_unrod_cnr = Universe(name='Rfl Assembly Unrodded Corner')
assembly_uo2_rod = Universe(name='UO2 Assembly Unrodded')
assembly_mox_rod = Universe(name='MOX Assembly Unrodded')
assembly_rfl_rod = Universe(name='Rfl Assembly Rodded')
assembly_rfl_unrod = Universe(name='Rfl Assembly Unrodded')

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
root_cell = Cell(name='Full Geometry')
root_cell.addSurface(halfspace=+1, surface=xmin)
root_cell.addSurface(halfspace=-1, surface=xmax)
root_cell.addSurface(halfspace=+1, surface=ymin)
root_cell.addSurface(halfspace=-1, surface=ymax)
root_cell.addSurface(halfspace=+1, surface=zmin)
root_cell.addSurface(halfspace=-1, surface=zmax)

root_universe = Universe(name='Root Universe')
root_universe.addCell(root_cell)
