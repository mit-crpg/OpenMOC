from openmoc import *
from openmoc.cuda import GPUSolver
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')

log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
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

# Create Circles for the fuel as well as to discretize the moderator into rings
fuel_radius = Circle(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = Circle(x=0.0, y=0.0, radius=0.62)
moderator_outer_radius = Circle(x=0.0, y=0.0, radius=0.58)


###############################################################################
######################   Creating Cells and Universes   #######################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator_ring1 = CellBasic(sectors=8)
moderator_ring2 = CellBasic(sectors=8)
moderator_ring3 = CellBasic(sectors=8)
moderator_ring1.setMaterial(materials['Water'])
moderator_ring2.setMaterial(materials['Water'])
moderator_ring3.setMaterial(materials['Water'])
moderator_ring1.addSurface(+1, fuel_radius)
moderator_ring1.addSurface(-1, moderator_inner_radius)
moderator_ring2.addSurface(+1, moderator_inner_radius)
moderator_ring2.addSurface(-1, moderator_outer_radius)
moderator_ring3.addSurface(+1, moderator_outer_radius)

# UO2 pin cell
uo2_cell = CellBasic(rings=3, sectors=8)
uo2_cell.setMaterial(materials['UO2'])
uo2_cell.addSurface(-1, fuel_radius)

uo2 = Universe(name='UO2')
uo2.addCell(uo2_cell)
uo2.addCell(moderator_ring1)
uo2.addCell(moderator_ring2)
uo2.addCell(moderator_ring3)

# 4.3% MOX pin cell
mox43_cell = CellBasic(rings=3, sectors=8)
mox43_cell.setMaterial(materials['MOX-4.3%'])
mox43_cell.addSurface(-1, fuel_radius)

mox43 = Universe(name='MOX-4.3%')
mox43.addCell(mox43_cell)
mox43.addCell(moderator_ring1)
mox43.addCell(moderator_ring2)
mox43.addCell(moderator_ring3)

# 7% MOX pin cell
mox7_cell = CellBasic(rings=3, sectors=8)
mox7_cell.setMaterial(materials['MOX-7%'])
mox7_cell.addSurface(-1, fuel_radius)

mox7 = Universe(name='MOX-7%')
mox7.addCell(mox7_cell)
mox7.addCell(moderator_ring1)
mox7.addCell(moderator_ring2)
mox7.addCell(moderator_ring3)

# 8.7% MOX pin cell
mox87_cell = CellBasic(rings=3, sectors=8)
mox87_cell.setMaterial(materials['MOX-8.7%'])
mox87_cell.addSurface(-1, fuel_radius)

mox87 = Universe(name='MOX-8.7%')
mox87.addCell(mox87_cell)
mox87.addCell(moderator_ring1)
mox87.addCell(moderator_ring2)
mox87.addCell(moderator_ring3)

# Fission chamber pin cell
fission_chamber_cell = CellBasic(rings=3, sectors=8)
fission_chamber_cell.setMaterial(materials['Fission Chamber'])
fission_chamber_cell.addSurface(-1, fuel_radius)

fission_chamber = Universe(name='Fission Chamber')
fission_chamber.addCell(fission_chamber_cell)
fission_chamber.addCell(moderator_ring1)
fission_chamber.addCell(moderator_ring2)
fission_chamber.addCell(moderator_ring3)

# Guide tube pin cell
guide_tube_cell = CellBasic(rings=3, sectors=8)
guide_tube_cell.setMaterial(materials['Guide Tube'])
guide_tube_cell.addSurface(-1, fuel_radius)

guide_tube = Universe(name='Guide Tube')
guide_tube.addCell(guide_tube_cell)
guide_tube.addCell(moderator_ring1)
guide_tube.addCell(moderator_ring2)
guide_tube.addCell(moderator_ring3)

# Reflector
reflector_cell = CellBasic(name='moderator')
reflector_cell.setMaterial(materials['Water'])

reflector = Universe(name='Reflector')
reflector.addCell(reflector_cell)

# CellFills
assembly1_cell = CellFill(name='Assembly 1')
assembly2_cell = CellFill(name='Assembly 2')
refined_reflector_cell = CellFill(name='Semi-Finely Spaced Reflector')
right_reflector_cell = CellFill(name='Right Reflector')
corner_reflector_cell = CellFill(name='Bottom Corner Reflector')
bottom_reflector_cell = CellFill(name='Bottom Reflector')

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
root_cell = CellFill(name='Full Geometry')
root_cell.addSurface(+1, boundaries[0])
root_cell.addSurface(-1, boundaries[1])
root_cell.addSurface(-1, boundaries[2])
root_cell.addSurface(+1, boundaries[3])

root_universe = Universe(name='Root Universe')
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
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
refined_reflector_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Right Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses(template)
right_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(name='Bottom Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bottom_reflector_cell.setFill(lattices[-1])

# 4 x 4 core to represent two bundles and water
lattices.append(Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)
lattices[-1].setUniverses([
     [assembly1,        assembly2,        right_reflector],
     [assembly2,        assembly1,        right_reflector],
     [bottom_reflector, bottom_reflector, corner_reflector]])
root_cell.setFill(lattices[-1])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = GPUSolver(geometry, track_generator)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
