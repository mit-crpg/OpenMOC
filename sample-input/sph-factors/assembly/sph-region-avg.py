import matplotlib.pyplot as plt
import numpy as np
import openmc.mgxs
import openmoc
from openmoc.materialize import compute_sph_factors, load_openmc_mgxs_lib
from openmoc.compatible import get_openmoc_geometry
import openmoc.plotter


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = openmoc.options.Options()

num_threads = options.getNumThreads()
spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                 Eigenvalue Calculation w/o SPH Factors
###############################################################################

# Initialize 2-group OpenMC multi-group cross section library for a pin cell
mgxs_lib = openmc.mgxs.Library.load_from_file(filename='region-avg')

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)

# Load cross section data
openmoc_materials = load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

# Initialize an OpenMOC TrackGenerator and Solver
openmoc_geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(openmoc_geometry, num_azim, spacing)
track_generator.generateTracks()

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue()
solver.printTimerReport()
keff_no_sph = solver.getKeff()


###############################################################################
#                Eigenvalue Calculation with SPH Factors
###############################################################################

# Compute SPH factors
sph, sph_mgxs_lib, = compute_sph_factors(mgxs_lib, max_fix_src_iters=50,
                                         max_domain_iters=1,
                                         track_spacing=spacing,
                                         num_azim=num_azim,
                                         num_threads=num_threads)

# Load the SPH-corrected MGXS library data
openmoc_materials = load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

# Initialize an OpenMOC TrackGenerator and Solver
openmoc_geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(openmoc_geometry, num_azim, spacing)
track_generator.generateTracks()

# Run an eigenvalue calculation with the SPH-corrected modifed MGXS library
solver.setTrackGenerator(track_generator)
solver.computeEigenvalue()
solver.printTimerReport()
keff_with_sph = solver.getKeff()

# Report the OpenMC and OpenMOC eigenvalues
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/o SPH: \t%1.5f', keff_no_sph)
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/ SPH: \t%1.5f', keff_with_sph)
openmoc.log.py_printf('RESULT', 'OpenMC keff: \t\t1.02555 +/- 0.00031')


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

# Plot the geometry
openmoc.plotter.plot_materials(openmoc_geometry)
openmoc.plotter.plot_cells(openmoc_geometry)
openmoc.plotter.plot_flat_source_regions(openmoc_geometry)

# Plot the scalar fluxes
energy_groups = list(range(1, mgxs_lib.num_groups+1))
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=energy_groups)
openmoc.plotter.plot_fission_rates(solver, transparent_zeros=True)

# Plot histograms of the SPH factors in each group
for group in range(sph_mgxs_lib.num_groups):
    fig = plt.figure()
    indices = np.where(sph[:,group] != 1.)
    plt.hist(sph[indices,group].flatten(), 20, normed=1, facecolor='green')
    plt.xlabel('SPH Factors')
    plt.ylabel('Probability')
    plt.title('Histogram of SPH Factors (Group {})'.format(group+1))
    plt.grid()
    plt.savefig('plots/sph-group-{}.png'.format(group+1), bbox_inches='tight')

# Create a scatter plot of the SPH factors in groups 1 and 2
fig = plt.figure()
indices = np.where(sph[:,0] != 1.)
plt.scatter(sph[indices,0].flatten(), sph[indices,1].flatten())
plt.xlabel('Group 1')
plt.xlabel('Group 2')
plt.title('SPH Factors')
plt.grid()
plt.savefig('plots/sph-scatter.png', bbox_inches='tight')
