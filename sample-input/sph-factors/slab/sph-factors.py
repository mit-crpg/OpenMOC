import openmoc
import openmc.mgxs

import numpy as np
import matplotlib

# Enable Matplotib to work for headless nodes
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                 Eigenvalue Calculation w/o SPH Factors
###############################################################################

# Initialize 2-group OpenMC multi-group cross section library for a pin cell
mgxs_lib = openmc.mgxs.Library.load_from_file(filename='mgxs', directory='.')

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = \
    openmoc.opencg_compatible.get_openmoc_geometry(mgxs_lib.opencg_geometry)

# Load cross section data
openmoc_materials = \
    openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

# Initialize an OpenMOC TrackGenerator and Solver
track_generator = openmoc.TrackGenerator(openmoc_geometry, opts.num_azim,
                                         opts.track_spacing)
track_generator.generateTracks()

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue()
solver.printTimerReport()
keff_no_sph = solver.getKeff()

# Extract the OpenMOC scalar fluxes
fluxes_no_sph = openmoc.process.get_scalar_fluxes(solver)


###############################################################################
#                Eigenvalue Calculation with SPH Factors
###############################################################################

# Compute SPH factors
sph, sph_mgxs_lib, sph_indices = \
    openmoc.materialize.compute_sph_factors(mgxs_lib, 
                                            track_spacing=opts.track_spacing,
                                            num_azim=opts.num_azim,
                                            num_threads=opts.num_omp_threads)

# Load the SPH-corrected MGXS library data
materials = \
    openmoc.materialize.load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

# Run an eigenvalue calculation with the SPH-corrected modifed MGXS library
solver.computeEigenvalue()
solver.printTimerReport()
keff_with_sph = solver.getKeff()

# Report the OpenMC and OpenMOC eigenvalues
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/o SPH: \t%1.5f', keff_no_sph)
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/ SPH: \t%1.5f', keff_with_sph)
openmoc.log.py_printf('RESULT', 'OpenMC keff: \t\t0.96281 +/- 0.00009')

# Extract the OpenMOC scalar fluxes
fluxes_sph = openmoc.process.get_scalar_fluxes(solver) * sph


###############################################################################
#                       Plottting Scalar Fluxes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

# Allocate arrays for SR-specific data to extract from OpenMOC model
num_srs = openmoc_geometry.getNumSRs()
cell_ids = np.zeros(num_srs, dtype=np.int)
centroids = np.zeros(num_srs, dtype=np.float)
volumes = np.zeros(num_srs, dtype=np.float)

# Find the cell IDs, volumes, centroids and fluxes for each SR
for sr_id in range(num_srs):
    cell = openmoc_geometry.findCellContainingSR(sr_id)
    cell_ids[sr_id] = cell.getId()
    volumes[sr_id] = solver.getSRVolume(sr_id)
    centroids[sr_id] = cell.getMinX()

# Organize cell IDs, volumes and fluxes in order of increasing centroid
indices = np.argsort(centroids)
centroids = centroids[indices]
cell_ids = cell_ids[indices]
volumes = volumes[indices]
fluxes_no_sph = fluxes_no_sph[indices,:]
fluxes_sph = fluxes_sph[indices,:]

# Get OpenMC fluxes
tot_fiss_src = 0.
openmc_fluxes = np.zeros((num_srs, mgxs_lib.num_groups))
for sr_id, cell_id in enumerate(cell_ids):

    # Get NuFissionXS for cell from MGXS Library
    mgxs = mgxs_lib.get_mgxs(cell_id, 'nu-fission')

    # Store this cell's flux
    flux_tally = mgxs.tallies['flux']
    openmc_fluxes[sr_id, :] = flux_tally.get_values().flatten()

    # Increment the total fission source
    nu_fission = mgxs.tallies['nu-fission']
    tot_fiss_src += np.sum(nu_fission.mean)

# Normalize OpenMC flux to total fission source * volume to compare to OpenMOC
openmc_fluxes /= volumes[:,np.newaxis] * tot_fiss_src

# Plot the OpenMOC and OpenMC spatially-varying fluxes
for group in range(mgxs_lib.num_groups):
    fig = plt.figure()
    plt.plot(centroids, openmc_fluxes[:,group])
    plt.plot(centroids, fluxes_no_sph[:,group])
    plt.plot(centroids, fluxes_sph[:,group])
    plt.legend(['openmc', 'openmoc (w/o sph)', 'openmoc (sph)'],loc='best')
    plt.title('Volume-Averaged Scalar Flux (Group {})'.format(group+1))
    plt.xlabel('x [cm]')
    plt.ylabel('flux')
    plt.savefig('flux-group-{}.png'.format(group+1), bbox_inches='tight')
