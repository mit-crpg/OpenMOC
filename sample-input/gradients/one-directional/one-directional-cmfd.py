import sys
sys.path.append('..')
from cube import geometry, root_cell, left, right, top, bottom
import openmoc

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')
openmoc.log.py_printf('TITLE', \
  'Simulating a one group homogeneous, one directional gradient...')

###############################################################################
###########################   Creating CMFD Mesh    ###########################
###############################################################################
openmoc.log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(25,1)

###############################################################################
#########################   Load the Cubic Geometry   #########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing cubic geometry...')

left.setBoundaryType(openmoc.VACUUM)
right.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.REFLECTIVE)

geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.azim_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(opts.tolerance)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()

###############################################################################
############################    Generating Plots   ############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=100)
openmoc.plotter.plot_cells(geometry, gridsize=100)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=100)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1])

openmoc.log.py_printf('TITLE', 'Finished')
