import openmoc
import openmc.openmoc_compatible
import openmc.mgxs

import numpy as np
import matplotlib

# Enable Matplotib to work for headless nodes
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


opts = openmoc.options.Options()
openmoc.log.set_log_level('NORMAL')


###############################################################################
#                 Eigenvalue Calculation w/o SPH Factors
###############################################################################

# Initialize 2-group OpenMC multi-group cross section library for a pin cell
mgxs_lib = openmc.mgxs.Library.load_from_file(filename='mgxs', directory='.')

# Create an OpenMOC Geometry from the OpenMOC Geometry
openmoc_geometry = \
    openmc.openmoc_compatible.get_openmoc_geometry(mgxs_lib.geometry)

# Load cross section data
openmoc_materials = \
    openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

# Initialize FSRs
openmoc_geometry.initializeFlatSourceRegions()

# Initialize an OpenMOC TrackGenerator
track_generator = openmoc.TrackGenerator(
    openmoc_geometry, opts.num_azim, opts.azim_spacing)
track_generator.generateTracks()

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()
keff_no_sph = solver.getKeff()

# Extract the OpenMOC scalar fluxes
fluxes_no_sph = openmoc.process.get_scalar_fluxes(solver)


###############################################################################
#                Eigenvalue Calculation with SPH Factors
###############################################################################

# Compute SPH factors
sph, sph_mgxs_lib, sph_indices = \
    openmoc.materialize.compute_sph_factors(
        mgxs_lib, azim_spacing=opts.azim_spacing,
        num_azim=opts.num_azim, num_threads=opts.num_omp_threads)

# Load the SPH-corrected MGXS library data
materials = \
    openmoc.materialize.load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

# Run an eigenvalue calculation with the SPH-corrected modified MGXS library
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()
keff_with_sph = solver.getKeff()

# Report the OpenMC and OpenMOC eigenvalues
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/o SPH: \t%1.5f', keff_no_sph)
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/ SPH: \t%1.5f', keff_with_sph)
openmoc.log.py_printf('RESULT', 'OpenMC keff: \t\t1.17574 +/- 0.00086')


###############################################################################
#                         Extracting Scalar Fluxes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

# Plot the cells
openmoc.plotter.plot_cells(openmoc_geometry)

# Extract the OpenMOC scalar fluxes
fluxes_sph = openmoc.process.get_scalar_fluxes(solver)
fluxes_sph *= sph

# Extract the OpenMC scalar fluxes
num_fsrs = openmoc_geometry.getNumFSRs()
num_groups = openmoc_geometry.getNumEnergyGroups()
openmc_fluxes = np.zeros((num_fsrs, num_groups), dtype=np.float64)
nufission_xs = np.zeros((num_fsrs, num_groups), dtype=np.float64)

# Get the OpenMC flux in each FSR
for fsr in range(num_fsrs):

    # Find the OpenMOC cell and volume for this FSR
    openmoc_cell = openmoc_geometry.findCellContainingFSR(fsr)
    cell_id = openmoc_cell.getId()
    fsr_volume = track_generator.getFSRVolume(fsr)

    # Store the volume-averaged flux
    mgxs = mgxs_lib.get_mgxs(cell_id, 'nu-fission')
    flux = mgxs.tallies['flux'].mean.flatten()
    flux = np.flipud(flux) / fsr_volume
    openmc_fluxes[fsr, :] = flux
    nufission_xs[fsr, :] = mgxs.get_xs(nuclide='all')

# Extract energy group edges
group_edges = mgxs_lib.energy_groups.group_edges
group_edges += 1e-3     # Adjust lower bound to 1e-3 eV (for loglog scaling)

# Compute difference in energy bounds for each group
group_edges = np.flipud(group_edges)

# Normalize fluxes with the fission source
openmc_fluxes /= np.sum(openmc_fluxes * nufission_xs)
fluxes_sph /= np.sum(fluxes_sph * nufission_xs)
fluxes_no_sph /= np.sum(fluxes_no_sph * nufission_xs)


###############################################################################
#                 Plot the OpenMC, OpenMOC Scalar Fluxes
###############################################################################

# Extend the mgxs values array for matplotlib's step plot of fluxes
openmc_fluxes = np.insert(openmc_fluxes, 0, openmc_fluxes[:,0], axis=1)
fluxes_no_sph = np.insert(fluxes_no_sph, 0, fluxes_no_sph[:,0], axis=1)
fluxes_sph = np.insert(fluxes_sph, 0, fluxes_sph[:,0], axis=1)

# Plot OpenMOC and OpenMC fluxes in each FSR
for fsr in range(num_fsrs):

    # Get the OpenMOC cell and material for this FSR
    cell = openmoc_geometry.findCellContainingFSR(fsr)
    material_name = cell.getFillMaterial().getName()

    # Create a step plot for the MGXS
    fig = plt.figure()
    plt.plot(group_edges, openmc_fluxes[fsr,:],
             drawstyle='steps', color='r', linewidth=2)
    plt.plot(group_edges, fluxes_no_sph[fsr,:],
             drawstyle='steps', color='b', linewidth=2)
    plt.plot(group_edges, fluxes_sph[fsr,:],
             drawstyle='steps', color='g', linewidth=2)

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux')
    plt.title('Normalized Flux ({0})'.format(material_name))
    plt.xlim((min(group_edges), max(group_edges)))
    plt.legend(['openmc', 'openmoc w/o sph', 'openmoc w/ sph'], loc='best')
    plt.grid()
    filename = 'plots/flux-{0}.png'.format(material_name.replace(' ', '-'))
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


###############################################################################
#                 Plot OpenMC-to-OpenMOC Scalar Flux Errors
###############################################################################

# Compute the percent relative error in the flux
rel_err_no_sph = np.zeros(openmc_fluxes.shape)
rel_err_sph = np.zeros(openmc_fluxes.shape)

for fsr in range(num_fsrs):
    delta_flux_no_sph = fluxes_no_sph[fsr,:] - openmc_fluxes[fsr,:]
    delta_flux_sph = fluxes_sph[fsr,:] - openmc_fluxes[fsr,:]
    rel_err_no_sph[fsr,:] = delta_flux_no_sph / openmc_fluxes[fsr,:] * 100.
    rel_err_sph[fsr,:] = delta_flux_sph / openmc_fluxes[fsr,:] * 100.

# Plot OpenMOC relative flux errors in each FSR
for fsr in range(num_fsrs):

    # Get the OpenMOC cell and material for this FSR
    cell = openmoc_geometry.findCellContainingFSR(fsr)
    material_name = cell.getFillMaterial().getName()


    # Create a step plot for the MGXS
    fig = plt.figure()
    plt.plot(group_edges, rel_err_no_sph[fsr,:],
             drawstyle='steps', color='r', linewidth=2)
    plt.plot(group_edges, rel_err_sph[fsr,:],
             drawstyle='steps', color='b', linewidth=2)

    plt.xscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Relative Error [%]')
    plt.title('OpenMOC-to-OpenMC Flux Rel. Err. ({0})'.format(material_name))
    plt.xlim((min(group_edges), max(group_edges)))
    plt.legend(['openmoc w/o sph', 'openmoc w/ sph'], loc='best')
    plt.grid()
    filename = 'plots/rel-err-{0}.png'.format(material_name.replace(' ', '-'))
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
