import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')

openmoc.log.py_printf('TITLE', \
  'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
#                              Creating Cmfd mesh
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,51)
#cmfd.setGroupStructure([1,4,8])
cmfd.setKNearest(3)

from geometry_ls import geometry
geometry.setCmfd(cmfd)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.track_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = openmoc.CPULSSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

#openmoc.process.compute_fission_rates(solver)

#openmoc.plotter.plot_materials(geometry, gridsize=250)
#openmoc.plotter.plot_cells(geometry, gridsize=250)
#openmoc.plotter.plot_cmfd_cells(geometry, cmfd, gridsize=250)
#openmoc.plotter.plot_flat_source_regions(geometry, gridsize=250)
#openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
#openmoc.plotter.plot_spatial_fluxes_ls(solver, energy_groups=[1,2,3,4,5,6,7])
#openmoc.plotter.plot_fission_rates(solver, gridsize=250, norm=True)

openmoc.log.py_printf('TITLE', 'Finished')