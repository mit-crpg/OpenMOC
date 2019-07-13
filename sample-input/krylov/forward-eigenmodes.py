import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()
num_modes = 5

openmoc.log.set_log_level('NORMAL')
openmoc.log.py_printf('TITLE', 'Computing %d forward eigenmodes', num_modes)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

from geometry import geometry
geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.azim_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

# Initialize a CPUSolver to perform forward fixed source calculations
cpu_solver = openmoc.CPUSolver(track_generator)
cpu_solver.setNumThreads(opts.num_omp_threads)

# Initialize IRAMSolver to perform forward eigenmode calculation
iram_solver = openmoc.krylov.IRAMSolver(cpu_solver)
iram_solver.computeEigenmodes(num_modes=num_modes, solver_mode=openmoc.FORWARD)

# Report the forward eigenvalues to the user
eigenvalues = iram_solver._eigenvalues
openmoc.log.py_printf('RESULT', 'Forward eigenvalues: %s', str(eigenvalues))


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=500)
openmoc.plotter.plot_cells(geometry, gridsize=500)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=500)
openmoc.plotter.plot_eigenmode_fluxes(iram_solver, gridsize=250,
                                      energy_groups=[1,2,3,4,5,6,7],
                                      eigenmodes=range(1,num_modes+1))

openmoc.log.py_printf('TITLE', 'Finished')
