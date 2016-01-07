from openmoc import *
import openmoc.log as log
import openmoc.materialize as materialize
from surfaces import *

log.set_log_level('NORMAL')
log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../../c5g7-materials.py')


###############################################################################
######################   Creating Cells and Universes   #######################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator = Cell()
moderator.setFill(materials['Water'])
moderator.addSurface(+1, fuel_radius)
moderator.setNumRings(2)
moderator.setNumSectors(8)

# UO2 pin cell
uo2_cell = Cell()
uo2_cell.setNumRings(3)
uo2_cell.setNumSectors(8)
uo2_cell.setFill(materials['UO2'])
uo2_cell.addSurface(-1, fuel_radius)

uo2 = Universe(name='UO2')
uo2.addCell(uo2_cell)
uo2.addCell(moderator)

# 4.3% MOX pin cell
mox43_cell = Cell()
mox43_cell.setNumRings(3)
mox43_cell.setNumSectors(8)
mox43_cell.setFill(materials['MOX-4.3%'])
mox43_cell.addSurface(-1, fuel_radius)

mox43 = Universe(name='MOX-4.3%')
mox43.addCell(mox43_cell)
mox43.addCell(moderator)

# 7% MOX pin cell
mox7_cell = Cell()
mox7_cell.setNumRings(3)
mox7_cell.setNumSectors(8)
mox7_cell.setFill(materials['MOX-7%'])
mox7_cell.addSurface(-1, fuel_radius)

mox7 = Universe(name='MOX-7%')
mox7.addCell(mox7_cell)
mox7.addCell(moderator)

# 8.7% MOX pin cell
mox87_cell = Cell()
mox87_cell.setNumRings(3)
mox87_cell.setNumSectors(8)
mox87_cell.setFill(materials['MOX-8.7%'])
mox87_cell.addSurface(-1, fuel_radius)

mox87 = Universe(name='MOX-8.7%')
mox87.addCell(mox87_cell)
mox87.addCell(moderator)

# Fission chamber pin cell
fission_chamber_cell = Cell()
fission_chamber_cell.setNumRings(3)
fission_chamber_cell.setNumSectors(8)
fission_chamber_cell.setFill(materials['Fission Chamber'])
fission_chamber_cell.addSurface(-1, fuel_radius)

fission_chamber = Universe(name='Fission Chamber')
fission_chamber.addCell(fission_chamber_cell)
fission_chamber.addCell(moderator)

# Guide tube pin cell
guide_tube_cell = Cell()
guide_tube_cell.setNumRings(3)
guide_tube_cell.setNumSectors(8)
guide_tube_cell.setFill(materials['Guide Tube'])
guide_tube_cell.addSurface(-1, fuel_radius)

guide_tube = Universe(name='Guide Tube')
guide_tube.addCell(guide_tube_cell)
guide_tube.addCell(moderator)

# Control rod pin cell
control_rod_cell = Cell()
control_rod_cell.setNumRings(3)
control_rod_cell.setNumSectors(8)
control_rod_cell.setFill(materials['Control Rod'])
control_rod_cell.addSurface(-1, fuel_radius)

control_rod = Universe(name='Control Rod')
control_rod.addCell(control_rod_cell)
control_rod.addCell(moderator)

# Reflector
reflector_cell = Cell(name='moderator')
reflector_cell.setFill(materials['Water'])

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
