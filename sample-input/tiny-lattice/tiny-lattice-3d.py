from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options
from geometry import *

###############################################################################
#                           Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.getNumThreads()
azim_spacing = options.getAzimSpacing()
num_azim = options.getNumAzimAngles()
polar_spacing = options.getPolarSpacing()
num_polar = options.getNumPolarAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, num_polar, azim_spacing, \
                                 polar_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

#plotter.plot_tracks_3d(track_generator)
#plotter.plot_segments_3d(track_generator)

log.py_printf('TITLE', 'Finished')
