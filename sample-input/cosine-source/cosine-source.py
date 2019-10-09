import openmoc
import math

use_gpu = False

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

boundary = openmoc.RectangularPrism(40., 40.)
boundary.setBoundaryType(openmoc.VACUUM)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

water_cell = openmoc.Cell(name='water')
water_cell.setFill(materials['Water'])

root_cell = openmoc.Cell(name='root cell')
root_cell.setRegion(boundary)


###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

water_univ = openmoc.Universe(name='water')
root_universe = openmoc.Universe(name='root universe')

water_univ.addCell(water_cell)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

# Number of lattice cells
num_x = 200
num_y = 200

# Compute widths of each lattice cell
width_x = (root_universe.getMaxX() - root_universe.getMinX()) / num_y
width_y = (root_universe.getMaxY() - root_universe.getMinY()) / num_x

# Create 2D array of Universes in each lattice cell
universes = [[[water_univ]*num_x for _ in range(num_y)]]

openmoc.log.py_printf('NORMAL', \
  'Creating a {0}x{0} lattice...'.format(num_x, num_y))

lattice = openmoc.Lattice(name='{0}x{1} lattice'.format(num_x, num_y))
lattice.setWidth(width_x=width_x, width_y=width_y)
lattice.setUniverses(universes)
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

# initialize the solver
if use_gpu:
    from openmoc.cuda import GPUSolver
    solver= GPUSolver(track_generator)
else:
    solver = openmoc.CPUSolver(track_generator)
    solver.setNumThreads(opts.num_omp_threads)
solver.initializeSolver(0)
solver.setConvergenceThreshold(opts.tolerance)

# Set the source in every cell to a cosine distribution
for fsr_id in range(solver.getGeometry().getNumFSRs()):

  # Get the coordinates of some point within the FSR
  pt = solver.getGeometry().getFSRPoint(fsr_id)
  x_pt = pt.getX()
  y_pt = pt.getY()

  # Set the FSR source for every group
  L = num_x * width_x / 2
  H = num_y * width_y / 2
  for g in range(materials['Water'].getNumEnergyGroups()):
    group = g + 1
    source_value = math.cos(x_pt/L) * math.cos(y_pt/H)
    solver.setFixedSourceByFSR(fsr_id, group, source_value)

  # NOTE: A more precise definition of the source would calculate the same
  # source values for all points within each flat source region. In this
  # example that is not the case. However, since the FSR discretization is
  # reasonably fine in this case, the slight error introduced from defining the
  # source based on one point in the FSR does not have a large impact on the
  # resulting flux shapes.

# Run the solver for the provided source
solver.computeFlux(opts.max_iters)
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
