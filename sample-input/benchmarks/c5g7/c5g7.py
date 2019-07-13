from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.num_omp_threads
azim_spacing = options.azim_spacing
num_azim = options.num_azim
tolerance = options.tolerance
max_iters = options.max_iters

log.set_log_level('NORMAL')

log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../../')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-32.13, name='left')
right = XPlane(x=32.13, name='right')
top = YPlane(y=32.13, name='top')
bottom = YPlane(y=-32.13, name='bottom')
left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(VACUUM)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(VACUUM)
boundaries = [left, right, top, bottom]

# Create ZCylinders for the fuel as well as to discretize the moderator into rings
fuel_radius = ZCylinder(x=0.0, y=0.0, radius=0.54)


###############################################################################
#                        Creating Cells and Universes
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator = Cell()
moderator.setNumSectors(8)
moderator.setNumRings(3)
moderator.setFill(materials['Water'])
moderator.addSurface(+1, fuel_radius)

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

# Reflector
reflector_cell = Cell(name='moderator')
reflector_cell.setFill(materials['Water'])

reflector = Universe(name='Reflector')
reflector.addCell(reflector_cell)

# Cells
assembly1_cell = Cell(name='Assembly 1')
assembly2_cell = Cell(name='Assembly 2')
refined_reflector_cell = Cell(name='Semi-Finely Spaced Reflector')
right_reflector_cell = Cell(name='Right Reflector')
corner_reflector_cell = Cell(name='Bottom Corner Reflector')
bottom_reflector_cell = Cell(name='Bottom Reflector')

assembly1 = Universe(name='Assembly 1')
assembly2 = Universe(name='Assembly 2')
refined_reflector = Universe(name='Semi-Finely Spaced Moderator')
right_reflector = Universe(name='Right Reflector')
corner_reflector = Universe(name='Bottom Corner Reflector')
bottom_reflector = Universe(name='Bottom Reflector')

assembly1.addCell(assembly1_cell)
assembly2.addCell(assembly2_cell)
refined_reflector.addCell(refined_reflector_cell)
right_reflector.addCell(right_reflector_cell)
corner_reflector.addCell(corner_reflector_cell)
bottom_reflector.addCell(bottom_reflector_cell)

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
lattices[-1].setUniverses([template])
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
lattices[-1].setUniverses([template])
assembly2_cell.setFill(lattices[-1])

# Sliced up water cells - semi finely spaced
lattices.append(Lattice(name='Semi-Finely Spaced Reflector'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
template = [[reflector] * 10] * 10
lattices[-1].setUniverses([template])
refined_reflector_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Right Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses([template])
right_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(name='Bottom Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses([template])
bottom_reflector_cell.setFill(lattices[-1])

# 4 x 4 core to represent two bundles and water
lattices.append(Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)
lattices[-1].setUniverses([[
     [assembly1,        assembly2,        right_reflector],
     [assembly2,        assembly1,        right_reflector],
     [bottom_reflector, bottom_reflector, corner_reflector]]])
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

track_generator = TrackGenerator(geometry, num_azim, azim_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=500)
plotter.plot_cells(geometry, gridsize=500)
plotter.plot_flat_source_regions(geometry, gridsize=500)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
