from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-5.0, name='left')
right = XPlane(x=5.0, name='right')
bottom = YPlane(y=-5.0, name='top')
top = YPlane(y=5.0, name='bottom')
boundaries = [left, right, top, bottom]

for boundary in boundaries: boundary.setBoundaryType(VACUUM)


###############################################################################
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

water_cell = CellBasic(name='water')
water_cell.setMaterial(materials['Water'])

fuel_cell = CellBasic(name='fuel')
fuel_cell.setMaterial(materials['UO2'])

source_cell = CellBasic(name='source')
source_cell.setMaterial(materials['Water'])

root_cell = CellFill(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)


###############################################################################
#                             Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

water_univ = Universe(name='water')
fuel_univ = Universe(name='fuel')
source_univ = Universe(name='source')
root_universe = Universe(name='root universe')

water_univ.addCell(water_cell)
fuel_univ.addCell(fuel_cell)
source_univ.addCell(source_cell)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

num_x = 100
num_y = 100
width_x = (root_universe.getMaxX() - root_universe.getMinX()) / num_y
width_y = (root_universe.getMaxY() - root_universe.getMinY()) / num_x
universes = [[water_univ]*num_x for _ in range(num_y)]

universes[15][15] = source_univ

for i in range(50,56):
  for j in range(50,56):
    universes[i][j] = fuel_univ

log.py_printf('NORMAL', 'Creating a {0}x{0} lattice...'.format(num_x, num_y))

lattice = Lattice(name='{0}x{1} lattice'.format(num_x, num_y))
lattice.setWidth(width_x=width_x, width_y=width_y)
lattice.setUniverses(universes)
root_cell.setFill(lattice)


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

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.setFixedSourceByCell(source_cell, group=3, source=1e-2)
solver.computeSource(max_iters, k_eff=0.6)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=250)
plotter.plot_cells(geometry, gridsize=250)
plotter.plot_flat_source_regions(geometry, gridsize=250)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
