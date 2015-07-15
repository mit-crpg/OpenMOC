import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = 2000

log.set_log_level('NORMAL')
log.py_printf('TITLE', 'Simulating a one group homogeneous, one directional'
    ' gradient...')

###############################################################################
#########################   Load the Cubic Geometry   #########################
###############################################################################

log.py_printf('NORMAL', 'Importing cubic geometry...')

import sys
sys.path.append('..')
import boundaries
boundaries.left_bc = VACUUM
boundaries.right_bc = REFLECTIVE
boundaries.top_bc = REFLECTIVE
boundaries.bottom_bc = REFLECTIVE
from cube import *

geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
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

###############################################################################
############################    Generating Plots   ############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=100)
plotter.plot_cells(geometry, gridsize=100)
plotter.plot_flat_source_regions(geometry, gridsize=100)
plotter.plot_spatial_fluxes(solver, energy_groups=[1])

log.py_printf('TITLE', 'Finished')
