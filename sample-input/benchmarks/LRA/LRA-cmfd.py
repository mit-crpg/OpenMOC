import openmoc
import openmoc.plotter as plotter
import openmoc.process as process
from openmoc.options import Options
from lattices import lattices, universes

###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

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

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(110, 110)
cmfd.setCentroidUpdateOn(False)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setCmfd(cmfd)
geometry.setRootUniverse(universes['Root'])


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

process.compute_fission_rates(solver)

openmoc.plotter.plot_materials(geometry, gridsize=500)
openmoc.plotter.plot_cells(geometry, gridsize=500)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=500)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2])

openmoc.log.py_printf('TITLE', 'Finished')
