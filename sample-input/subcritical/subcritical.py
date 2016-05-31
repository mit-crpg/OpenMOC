import numpy as np
import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

left =openmoc. XPlane(x=-5.0, name='left')
right = openmoc.XPlane(x=5.0, name='right')
bottom = openmoc.YPlane(y=-5.0, name='bottom')
top = openmoc.YPlane(y=5.0, name='top')
boundaries = [left, right, top, bottom]

for boundary in boundaries: boundary.setBoundaryType(openmoc.VACUUM)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

water_cell = openmoc.Cell(name='water')
water_cell.setFill(materials['Water'])

fuel_cell = openmoc.Cell(name='fuel')
fuel_cell.setFill(materials['UO2'])

source_cell = openmoc.Cell(name='source')
source_cell.setFill(materials['Water'])

root_cell = openmoc.Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)


###############################################################################
#                             Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

water_univ = openmoc.Universe(name='water')
fuel_univ = openmoc.Universe(name='fuel')
source_univ = openmoc.Universe(name='source')
root_universe = openmoc.Universe(name='root universe')

water_univ.addCell(water_cell)
fuel_univ.addCell(fuel_cell)
source_univ.addCell(source_cell)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

# Number of lattice cells
num_x = 100
num_y = 100

# Compute widths of each lattice cell
width_x = (root_universe.getMaxX() - root_universe.getMinX()) / num_y
width_y = (root_universe.getMaxY() - root_universe.getMinY()) / num_x

# Create 2D array of Universes in each lattice cell
universes = [[water_univ]*num_x for _ in range(num_y)]

# Place fixed source Universe at (x=10, y=10)
source_x = 2.5
source_y = 2.5
lat_x = (root_universe.getMaxX() - source_x) / width_x
lat_y = (root_universe.getMaxY() - source_y) / width_y
universes[int(lat_x)][int(lat_y)] = source_univ

# Place fuel Universes at (x=[-0.5, 0.5], y=[-0.5,0.5])
fuel_x = np.array([-0.5,0.5])
fuel_y = np.array([-0.5,0.5])
lat_x = (root_universe.getMaxX() - fuel_x) / width_x
lat_y = (root_universe.getMaxY() - fuel_y) / width_y

for i in np.arange(lat_x[1], lat_x[0], dtype=np.int):
  for j in np.arange(lat_y[1], lat_y[0], dtype=np.int):
    universes[i][j] = fuel_univ

openmoc.log.py_printf('NORMAL', \
    'Creating a {0}x{0} lattice...'.format(num_x, num_y))

lattice = openmoc.Lattice(name='{0}x{1} lattice'.format(num_x, num_y))
lattice.setWidth(width_x=width_x, width_y=width_y)
lattice.setUniverses([universes])
root_cell.setFill(lattice)


###############################################################################
#                         Creating the Geometry
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.azim_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(opts.tolerance)
solver.setFixedSourceByCell(source_cell, group=3, source=1e-2)
solver.computeSource(opts.max_iters, k_eff=0.6)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=250)
openmoc.plotter.plot_cells(geometry, gridsize=250)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=250)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

openmoc.log.py_printf('TITLE', 'Finished')
