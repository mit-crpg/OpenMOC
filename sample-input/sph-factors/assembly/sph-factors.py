import copy
import numpy as np

import openmc.mgxs
import openmoc
from openmoc.materialize import compute_sph_factors, load_openmc_mgxs_lib
from openmoc.compatible import get_openmoc_geometry
import openmoc.plotter
from infermc import differentiate_mgxs_lib


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
#                 Construct OpenMOC Model OpenMC MGXS Library
###############################################################################

# Initialize 2-group OpenMC multi-group cross section library for a pin cell
mgxs_lib = openmc.mgxs.Library.load_from_file()

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)

# Load the OpenMC MGXS library data from the OpenMC MGXS library
diff_mgxs_lib = mgxs_lib.get_subdomain_avg_library()
openmoc_materials = load_openmc_mgxs_lib(diff_mgxs_lib, openmoc_geometry)

# Initialize an OpenMOC TrackGenerator and Solver
openmoc_geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(openmoc_geometry, num_azim, spacing)
track_generator.generateTracks()

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)


###############################################################################
#                 Eigenvalue Calculation w/o SPH Factors
###############################################################################

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue()
solver.printTimerReport()
keff_no_sph = solver.getKeff()


###############################################################################
#     Differentiate Fuel Pins with OpenCG Local Neighbor Symmetry (LNS)
###############################################################################

opencg_geometry = copy.deepcopy(mgxs_lib.opencg_geometry)
opencg_geometry.build_neighbors()
opencg_geometry.count_neighbors(first_level=0)

regions_to_clusters = np.zeros(opencg_geometry.num_regions, dtype=np.int)
counter = 1
for region in range(opencg_geometry.num_regions):
    coords = opencg_geometry.find_region(region)

    if 'Fuel' in coords.tail_node.cell.fill.name:
#        regions_to_clusters[region] = counter
#        counter += 1
        regions_to_clusters[region] = \
            opencg_geometry.get_neighbors_hash(region)

clusters_to_regions = opencg_geometry.differentiate(regions_to_clusters)


###############################################################################
#                Eigenvalue Calculation with SPH Factors
###############################################################################

# Compute SPH factors
'''
sph, sph_mgxs_lib, = compute_sph_factors(mgxs_lib, num_azim=num_azim, 
                                         max_fix_src_iters=50,
                                         mode='fissionable',
                                         max_domain_iters=1,
                                         track_spacing=spacing,
                                         num_threads=num_threads)
'''

#sph_mgxs_lib = diff_mgxs_lib
sph_mgxs_lib = \
    differentiate_mgxs_lib(mgxs_lib, clusters_to_regions, opencg_geometry)


# Load the SPH-corrected MGXS library data
openmoc_materials = load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

# Run an eigenvalue calculation with the SPH-corrected modifed MGXS library
solver.computeEigenvalue()
solver.printTimerReport()
keff_with_sph = solver.getKeff()

# Report the OpenMC and OpenMOC eigenvalues
sp = openmc.StatePoint(mgxs_lib.sp_filename)
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/o SPH: \t%1.5f', keff_no_sph)
openmoc.log.py_printf('RESULT', 'OpenMOC keff w/ SPH: \t%1.5f', keff_with_sph)
openmoc.log.py_printf('RESULT', 'OpenMC keff: \t\t%1.5f +/- %1.5f',
                                sp.k_combined[0], sp.k_combined[1])


openmoc.plotter.plot_materials(openmoc_geometry)
openmoc.plotter.plot_cells(openmoc_geometry)
openmoc.plotter.plot_flat_source_regions(openmoc_geometry)

# Use OpenMOC to plot the scalar fluxes
energy_groups = list(range(1, mgxs_lib.num_groups+1))
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=energy_groups)
