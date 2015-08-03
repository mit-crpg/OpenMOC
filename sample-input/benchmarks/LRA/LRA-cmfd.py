import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = openmoc.options.Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

openmoc.log.set_log_level('NORMAL')

openmoc.log.py_printf('TITLE', 'Simulating the LRA Benchmark Problem...')


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating Cmfd mesh...')

from geometry import geometry
cmfd = openmoc.Cmfd()
cmfd.setLatticeStructure(110,110)
geometry.setCmfd(cmfd)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=500)
openmoc.plotter.plot_cells(geometry, gridsize=500)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=500)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2])
openmoc.plotter.plot_cmfd_cells(geometry, cmfd, gridsize=500)

openmoc.log.py_printf('TITLE', 'Finished')

