from openmoc import *
from openmoc.krylov import IRAMSolver
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options
import numpy as np


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()
num_modes = 5

log.set_log_level('NORMAL')

log.py_printf('TITLE', 'Computing %d eigenmodes', num_modes)


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-64.26, name='left')
right = XPlane(x=64.26, name='right')
top = YPlane(y=64.26, name='top')
bottom = YPlane(y=-64.26, name='bottom')
left.setBoundaryType(VACUUM)
right.setBoundaryType(VACUUM)
top.setBoundaryType(VACUUM)
bottom.setBoundaryType(VACUUM)
boundaries = [left, right, top, bottom]

# Create Circles for the fuel as well as to discretize the moderator into rings
fuel_radius = Circle(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = Circle(x=0.0, y=0.0, radius=0.62)
moderator_outer_radius = Circle(x=0.0, y=0.0, radius=0.58)


###############################################################################
#                        Creating Cells and Universes
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator_ring1 = Cell()
moderator_ring2 = Cell()
moderator_ring3 = Cell()
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
uo2_cell = Cell()
uo2_cell.setNumRings(3)
uo2_cell.setNumSectors(8)
uo2_cell.setFill(materials['UO2'])
uo2_cell.addSurface(-1, fuel_radius)

uo2 = Universe(name='UO2')
uo2.addCell(uo2_cell)
uo2.addCell(moderator_ring1)
uo2.addCell(moderator_ring2)
uo2.addCell(moderator_ring3)

# 4.3% MOX pin cell
mox43_cell = Cell()
mox43_cell.setNumRings(3)
mox43_cell.setNumSectors(8)
mox43_cell.setFill(materials['MOX-4.3%'])
mox43_cell.addSurface(-1, fuel_radius)

mox43 = Universe(name='MOX-4.3%')
mox43.addCell(mox43_cell)
mox43.addCell(moderator_ring1)
mox43.addCell(moderator_ring2)
mox43.addCell(moderator_ring3)

# 7% MOX pin cell
mox7_cell = Cell()
mox7_cell.setNumRings(3)
mox7_cell.setNumSectors(8)
mox7_cell.setFill(materials['MOX-7%'])
mox7_cell.addSurface(-1, fuel_radius)

mox7 = Universe(name='MOX-7%')
mox7.addCell(mox7_cell)
mox7.addCell(moderator_ring1)
mox7.addCell(moderator_ring2)
mox7.addCell(moderator_ring3)

# 8.7% MOX pin cell
mox87_cell = Cell()
mox87_cell.setNumRings(3)
mox87_cell.setNumSectors(8)
mox87_cell.setFill(materials['MOX-8.7%'])
mox87_cell.addSurface(-1, fuel_radius)

mox87 = Universe(name='MOX-8.7%')
mox87.addCell(mox87_cell)
mox87.addCell(moderator_ring1)
mox87.addCell(moderator_ring2)
mox87.addCell(moderator_ring3)

# Fission chamber pin cell
fission_chamber_cell = Cell()
fission_chamber_cell.setNumRings(3)
fission_chamber_cell.setNumSectors(8)
fission_chamber_cell.setFill(materials['Fission Chamber'])
fission_chamber_cell.addSurface(-1, fuel_radius)

fission_chamber = Universe(name='Fission Chamber')
fission_chamber.addCell(fission_chamber_cell)
fission_chamber.addCell(moderator_ring1)
fission_chamber.addCell(moderator_ring2)
fission_chamber.addCell(moderator_ring3)

# Guide tube pin cell
guide_tube_cell = Cell()
guide_tube_cell.setNumRings(3)
guide_tube_cell.setNumSectors(8)
guide_tube_cell.setFill(materials['Guide Tube'])
guide_tube_cell.addSurface(-1, fuel_radius)

guide_tube = Universe(name='Guide Tube')
guide_tube.addCell(guide_tube_cell)
guide_tube.addCell(moderator_ring1)
guide_tube.addCell(moderator_ring2)
guide_tube.addCell(moderator_ring3)

# Reflector
refl_cell = Cell(name='moderator')
refl_cell.setFill(materials['Water'])

reflector = Universe(name='Reflector')
reflector.addCell(refl_cell)

# Cells
assembly1_cell = Cell(name='Assembly 1')
assembly2_cell = Cell(name='Assembly 2')
refined_refl_cell = Cell(name='Semi-Finely Spaced Reflector')
right_refl_cell = Cell(name='Right Reflector')
left_refl_cell = Cell(name='Left Reflector')
bot_right_refl_cell = Cell(name='Bottom Right Corner Reflector')
top_right_refl_cell = Cell(name='Top Right Corner Reflector')
bot_left_refl_cell = Cell(name='Bottom Left Corner Reflector')
top_left_refl_cell = Cell(name='Top Left Corner Reflector')
bot_refl_cell = Cell(name='Bottom Reflector')
top_refl_cell = Cell(name='Top Reflector')

assembly1 = Universe(name='Assembly 1')
assembly2 = Universe(name='Assembly 2')
refined_refl = Universe(name='Semi-Finely Spaced Moderator')
right_refl = Universe(name='Right Reflector')
left_refl = Universe(name='Left Reflector')
bot_right_refl = Universe(name='Bottom Right Corner Reflector')
top_right_refl = Universe(name='Top Right Corner Reflector')
bot_left_refl = Universe(name='Bottom Left Corner Reflector')
top_left_refl = Universe(name='Top Left Corner Reflector')
bot_refl = Universe(name='Bottom Reflector')
top_refl = Universe(name='Top Reflector')

assembly1.addCell(assembly1_cell)
assembly2.addCell(assembly2_cell)
refined_refl.addCell(refined_refl_cell)
right_refl.addCell(right_refl_cell)
left_refl.addCell(left_refl_cell)
bot_right_refl.addCell(bot_right_refl_cell)
top_right_refl.addCell(top_right_refl_cell)
bot_left_refl.addCell(bot_left_refl_cell)
top_left_refl.addCell(top_left_refl_cell)
bot_refl.addCell(bot_refl_cell)
top_refl.addCell(top_refl_cell)

# Root Cell/Universe
root_cell = Cell(name='Full Geometry')
root_cell.addSurface(+1, boundaries[0])
root_cell.addSurface(-1, boundaries[1])
root_cell.addSurface(-1, boundaries[2])
root_cell.addSurface(+1, boundaries[3])

root_universe = Universe(name='Root Universe')
root_universe.addCell(root_cell)


###############################################################################
#                             Creating Lattices
###############################################################################

log.py_printf('NORMAL', 'Creating lattices...')

lattices = list()

# Top left, bottom right 17 x 17 assemblies
lattices.append(Lattice(name='Assembly 1'))
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
lattices[-1].setUniverses(template)
assembly1_cell.setFill(lattices[-1])

# Top right, bottom left 17 x 17 assemblies
lattices.append(Lattice(name='Assembly 2'))
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
lattices[-1].setUniverses(template)
assembly2_cell.setFill(lattices[-1])

# Sliced up water cells - semi finely spaced
lattices.append(Lattice(name='Semi-Finely Spaced Reflector'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
template = [[reflector] * 10] * 10
lattices[-1].setUniverses(template)
refined_refl_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Right Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_refl] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses(template)
right_refl_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Left Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 6 +  [refined_refl] * 11] * 17
lattices[-1].setUniverses(template)
left_refl_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Right Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_refl] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bot_right_refl_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Top Right Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[refined_refl] * 11 + [reflector] * 6] * 11
lattices[-1].setUniverses(template)
top_right_refl_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Left Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 6 + [refined_refl] * 11] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bot_left_refl_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Top Left Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[reflector] * 6 + [refined_refl] * 11] * 11
lattices[-1].setUniverses(template)
top_left_refl_cell.setFill(lattices[-1])

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(name='Bottom Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_refl] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bot_refl_cell.setFill(lattices[-1])

# Sliced up water cells for top of geometry
lattices.append(Lattice(name='Top Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[refined_refl] * 17] * 11
lattices[-1].setUniverses(template)
top_refl_cell.setFill(lattices[-1])

# 4 x 4 core to represent two bundles and water
lattices.append(Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)

lattices[-1].setUniverses([
  [top_left_refl, top_refl,  top_refl,  top_refl,  top_refl,  top_right_refl],
  [left_refl,     assembly1, assembly2, assembly2, assembly1, right_refl    ],
  [left_refl,     assembly2, assembly1, assembly1, assembly2, right_refl    ],
  [left_refl,     assembly2, assembly1, assembly1, assembly2, right_refl    ],
  [left_refl,     assembly1, assembly2, assembly2, assembly1, right_refl    ],
  [bot_left_refl, bot_refl,  bot_refl,  bot_refl,  bot_refl,  bot_right_refl]])
root_cell.setFill(lattices[-1])


###############################################################################
#                         Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

# Initialize a CPUSolver to perform fixed source calculations
cpu_solver = CPUSolver(track_generator)
cpu_solver.setNumThreads(num_threads)

# Initialize IRAMSolver to perform eigenmode calculation
iram_solver = IRAMSolver(cpu_solver)
iram_solver.computeEigenmodes(num_modes=num_modes)

# Report the eigenvalues to the user
eigenvalues = iram_solver._eigenvalues
log.py_printf('RESULT', 'The eigenvalues are: %s', str(eigenvalues))
        

###############################################################################
#                             Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=500)
plotter.plot_cells(geometry, gridsize=500)
plotter.plot_flat_source_regions(geometry, gridsize=500)
plotter.plot_eigenmode_fluxes(iram_solver, energy_groups=[1,2,3,4,5,6,7], gridsize=250)

log.py_printf('TITLE', 'Finished')
