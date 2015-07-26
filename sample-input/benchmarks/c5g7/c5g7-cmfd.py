import openmoc
import openmoc.log as log
import openmoc.plotter as plotter
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

log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
#                              Creating Cmfd mesh
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setMOCRelaxationFactor(0.6)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,51)
cmfd.setGroupStructure([1,4,8])

from geometry import geometry
geometry.setCmfd(cmfd)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

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

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=250)
plotter.plot_cells(geometry, gridsize=250)
plotter.plot_cmfd_cells(geometry, cmfd, gridsize=250)
plotter.plot_flat_source_regions(geometry, gridsize=250)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
plotter.plot_fission_rates(solver, gridsize=250)

log.py_printf('TITLE', 'Finished')
