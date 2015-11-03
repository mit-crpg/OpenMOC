##
# @file materialize.py
# @package openmoc.materialize
# @brief The materialize module provides utility functions to read and write
#        multi-group materials cross section data from HDF5 binary file format.
# @author William Boyd (wboyd@mit.edu)
# @date October 9, 2015

import sys
import os
import copy

import numpy as np

import openmoc

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf


def _get_domain(domains, domain_spec):
    """A helper routine to find materials/cells for load_from_hdf5(...)"""

    # If domain_spec is an integer, it must be a domain ID
    if isinstance(domain_spec, int) and domain_spec in domains:
        return domains[domain_spec]

    # If domain_spec is a string, it must be a domain name
    elif isinstance(domain_spec, str):
        for domain_id, domain in domains.items():
            if domain_spec == domain.getName():
                return domain

    # If domain could not be found
    return None


def _get_numpy_array(hdf5_group, key, suffix):
    """A helper routine to ensure that MGXS data is a proper NumPy array"""

    sigma = np.array(hdf5_group['{}/'.format(key) + suffix][...])
    sigma = np.atleast_1d(sigma)
    sigma = sigma.flatten()
    return sigma


##
# @brief This routine loads an HDF5 file of multi-group cross section data.
# @details The routine instantiates Material objects with multi-group cross
#          section data and returns a dictionary of each Material object keyed
#          by its ID. An OpenMOC Geometry may optionally be given and the
#          routine will directly insert the multi-group cross sections into each
#          Material in the Geometry. If a Geometry is passed in, Material
#          objects from the Geometry will be used in place of those instantiated
#          by this routine. An optional parameter for the domain types may be
#          used to define whether multi-group cross sections are tabulated by
#          material or cell in the HDF5 binary file.
# @param filename filename for cross sections HDF5 file (default is 'mgxs.h5')
# @param directory directory for cross sections HDF5 file (default is 'mgxs')
# @param geometry an optional geometry populated with materials, cells, etc.
# @param domain_type the domain type ('material' or 'cell') upon which the
#        cross sections are defined (default is 'material')
# @param suffix an optional string suffix to index the HDF5 file beyond the
#        assumed domain_type/domain_id/mgxs_type group sequence (default is '')
# @return a dictionary of Material objects keyed by ID
def load_from_hdf5(filename='mgxs.h5', directory='mgxs',
                   geometry=None, domain_type='material', suffix=''):

    # Create a h5py file handle for the file
    import h5py
    filename = os.path.join(directory, filename)
    f = h5py.File(filename, 'r')

    # Check that the file has an 'energy groups' attribute
    if '# groups' not in f.attrs:
        py_printf('ERROR', 'Unable to load HDF5 file "%s" since it does '
                           'not contain an \'# groups\' attribute', filename)

    if domain_type not in f.keys():
        py_printf('ERROR', 'Unable to load HDF5 file "%s" since it does '
                           'not contain domain type "%s"', filename, domain_type)

    if geometry and 'openmoc.Geometry' not in str(type(geometry)):
        py_printf('ERROR', 'Unable to load HDF5 file "%s" for "%s" which is not '
                           'an OpenMOC Geometry', filename, str(type(geometry)))

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    num_groups = int(f.attrs['# groups'])

    # If a Geometry was passed in, extract all cells or materials from it
    if geometry:
        if domain_type == 'material':
            domains = geometry.getAllMaterials()
        elif domain_type == 'cell':
            domains = geometry.getAllMaterialCells()
        else:
            py_printf('ERROR', 'Domain type "%s" is not supported', domain_type)

    # Iterate over all domains (e.g., materials or cells) in the HDF5 file
    for domain_spec in f[domain_type]:

        py_printf('INFO', 'Importing cross sections for %s "%s"',
                          domain_type, str(domain_spec))

        # Create shortcut to HDF5 group for this domain
        domain_group = f[domain_type][domain_spec]

        # If domain_spec is an integer, it is an ID; otherwise a string name
        if domain_spec.isdigit():
            domain_spec = int(domain_spec)
        else:
            domain_spec = str(domain_spec)

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:

            if domain_type == 'material':
                material = _get_domain(domains, domain_spec)

            elif domain_type == 'cell':
                cell = _get_domain(domains, domain_spec)
                material = cell.getFillMaterial()

                # If the user filled the Cell with a Material, clone it
                if material != None:
                    material = material.clone()

                # If the Cell does not contain a Material, create one for it
                else:
                    if isinstance(domain_spec, int):
                        material = openmoc.Material(id=domain_spec)
                    else:
                        material = openmoc.Material(name=domain_spec)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            if isinstance(domain_spec, int):
                material = openmoc.Material(id=domain_spec)
            else:
                material = openmoc.Material(name=domain_spec)

        # Add material to the collection
        materials[domain_spec] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in domain_group:
            sigma = _get_numpy_array(domain_group, 'transport', suffix)
            material.setSigmaT(sigma)
            py_printf('INFO', 'Loaded "transport" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        elif 'total' in domain_group:
            sigma = _get_numpy_array(domain_group, 'total', suffix)
            material.setSigmaT(sigma)
            py_printf('INFO', 'Loaded "total" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                                 '"%s %s"', domain_type, str(domain_spec))

        # Search for the fission production cross section
        if 'nu-fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-fission', suffix)
            material.setNuSigmaF(sigma)
            py_printf('INFO', 'Loaded "nu-fission" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                                 '"%s %s"', domain_type, str(domain_spec))

        # Search for the scattering matrix cross section
        if 'nu-scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-scatter matrix', suffix)
            material.setSigmaS(sigma)
            py_printf('INFO', 'Loaded "nu-scatter matrix" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        elif 'scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'scatter matrix', suffix)
            material.setSigmaS(sigma)
            py_printf('INFO', 'Loaded "scatter matrix" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "scatter matrix" found for "%s %s"', 
                      domain_type, str(domain_spec))

        # Search for chi (fission spectrum)
        if 'chi' in domain_group:
            chi = _get_numpy_array(domain_group, 'chi', suffix)
            material.setChi(chi)
            py_printf('INFO', 'Loaded "chi" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %s"',
                                 domain_type, str(domain_spec))

        # Search for optional cross sections
        if 'fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'fission', suffix)
            material.setSigmaF(sigma)
            py_printf('INFO', 'Loaded "fission" MGXS for "%s %s"', 
                      domain_type, str(domain_spec))

    # Return collection of materials
    return materials


##
# @brief This routine loads an OpenMC Library of multi-group cross section data.
# @details The routine instantiates Material objects with multi-group cross
#          section data and returns a dictionary of each Material object keyed
#          by its ID. An OpenMOC Geometry may optionally be given and the
#          routine will directly insert the multi-group cross sections into each
#          Material in the Geometry. If a Geometry is passed in, Material
#          objects from the Geometry will be used in place of those instantiated
#          by this routine.
# @param mgxs_lib an openmc.mgxs.Library object with cross section data
# @param geometry an optional geometry populated with materials, cells, etc.
# @return a dictionary of Material objects keyed by ID
def load_openmc_mgxs_lib(mgxs_lib, geometry=None):

    # Attempt to import openmc
    try:
        import openmc
    except ImportError:
        py_printf('ERROR', 'The OpenMC code must be installed on your system')

    if not isinstance(mgxs_lib, openmc.mgxs.Library):
        py_printf('ERROR', 'Unable to load cross sections from %s which is not '
                           'an openmc.mgxs.Library object', str(type(mgxs_lib)))

    if geometry and 'openmoc.Geometry' not in str(type(geometry)):
        py_printf('ERROR', 'Unable to load cross sections for "%s" which '
                           'is not an OpenMOC Geometry', str(type(geometry)))

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    num_groups = mgxs_lib.num_groups
    domain_type = mgxs_lib.domain_type

    # If a Geometry was passed in, extract all cells or materials from it
    if geometry:
        if domain_type == 'material':
            domains = geometry.getAllMaterials()
        elif domain_type == 'cell':
            domains = geometry.getAllMaterialCells()
        else:
            py_printf('ERROR', 'Unable to load a cross sections library with '
                               'domain type %s', mgxs_lib.domain_type)

    # Iterate over all domains (e.g., materials or cells) in the HDF5 file
    for domain in mgxs_lib.domains:

        py_printf('INFO', 'Importing cross sections for %s "%d"',
                          domain_type, domain.id)

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:

            if domain_type == 'material':
                material = _get_domain(domains, domain.id)

            elif domain_type == 'cell':
                cell = _get_domain(domains, domain.id)
                material = cell.getFillMaterial()

                # If the user filled the Cell with a Material, clone it
                if material != None:
                    # FIXME!!!
#                    material = material.clone()
                    material = material

                # If the Cell does not contain a Material, create one for it
                else:
                    material = openmoc.Material(id=domain)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            material = openmoc.Material(id=domain.id)

        # Add material to the collection
        materials[domain.id] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'transport')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
            py_printf('INFO', 'Loaded "transport" MGXS for "%s %d"', 
                      domain_type, domain.id)
        elif 'total' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'total')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
            py_printf('INFO', 'Loaded "total" MGXS for "%s %d"', 
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                                 '"%s %d"', domain_type, domain.id)

        # Search for the fission production cross section
        if 'nu-fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setNuSigmaF(sigma)
            py_printf('INFO', 'Loaded "nu-fission" MGXS for "%s %d"', 
                               domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                                 '"%s %d"', domain_type, domain.id)

        # Search for the scattering matrix cross section
        if 'nu-scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
            py_printf('INFO', 'Loaded "nu-scatter matrix" MGXS for "%s %d"', 
                      domain_type, domain.id)
        elif 'scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'scatter matrix').flatten()
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaS(sigma)
            py_printf('INFO', 'Loaded "scatter matrix" MGXS for "%s %d"', 
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "scatter matrix" or "nu-scatter matrix" '
                                 'found for "%s %d"', domain_type, domain.id)

        # Search for chi (fission spectrum)
        if 'chi' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'chi')
            chi = mgxs.get_xs(nuclides='sum')
            material.setChi(chi)
            py_printf('INFO', 'Loaded "chi" MGXS for "%s %d"', 
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %d"',
                                 domain_type, domain.id)

        # Search for optional cross sections
        if 'fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaF(sigma)
            py_printf('INFO', 'Loaded "fission" MGXS for "%s %d"', 
                      domain_type, domain.id)

    # Return collection of materials
    return materials


##
# @brief Compute SPH factors for an OpenMC multi-group cross section library
# @details This routine coputes SuPerHomogenisation (SPH) factors for an
#          an OpenMC MGXS library. The SPH scheme is outlined by Alain Hebert
#          in the following paper:
#
#          Hebert, A., "A Consistent Technique for the Pin-by-Pin
#          Homogenization of a Pressurized Water Reactor Assembly."
#          Nuclear Science and Engineering, 113 (3), pp. 227-233, 1993.
#
#          The SPH factors are needed to preserve reaction rates in
#          heterogeneous geometries. The energy condensation process leads
#          to a bias between ultrafine and coarse energy group calculations.
#          This bias is a result of the use of scalar flux-weighting to compute
#          MGXS without properly accounting for angular-dependence of the flux.
# @param mgxs_lib an openmc.mgxs.Library object
# @param max_iters the maximum number of SPH iterations
# @param sph_tol the tolerance on the SPH factor convergence
# @param fix_src_tol the tolerance on the MOC fixed source calculations
# @param num_azim number of azimuthal angles to use in fixed source calculations
# @param track_spacing the track spacing to use in fixed source calculations
# @return a NumPy array of SPH factors and a new openmc.mgxs.Library object
#         with the SPH factors applied to each MGXS object
def compute_sph_factors(mgxs_lib, max_iters=10, sph_tol=1E-6, fix_src_tol=1E-5,
                        num_azim=4, track_spacing=0.1, zcoord=0.0):

    import openmc.mgxs

    # For Python 2.X.X
    if sys.version_info[0] == 2:
        from compatible import get_openmoc_geometry
        from process import get_scalar_fluxes
    # For Python 3.X.X
    else:
        from openmoc.compatible import get_openmoc_geometry
        from openmoc.process import get_scalar_fluxes

    if not isinstance(mgxs_lib, openmc.mgxs.Library):
        py_printf('ERROR', 'Unable to compute SPH factors with %s which is '
                  'not an openmc.mgxs.Library', str(openmc.mgxs.Library))

    if mgxs_lib.domain_type not in ['material', 'cell']:
        py_printf('ERROR', 'SPH factors cannot be applied for an OpenMC MGXS '
                           'library of domain type %s', mgxs_lib.domain_type)

    # Create an OpenMOC Geometry from the OpenCG Geometry
    geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)

    # Load the MGXS library data into the OpenMOC geometry
    load_openmc_mgxs_lib(mgxs_lib, geometry)

    # Initialize an OpenMOC TrackGenerator
    geometry.initializeFlatSourceRegions()
    track_generator = openmoc.TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.setZCoord(zcoord)
    track_generator.generateTracks()

    # Initialize an OpenMOC Solver
    solver = openmoc.CPUSolver(track_generator)
    solver.setConvergenceThreshold(fix_src_tol)
    openmc_fluxes = _load_openmc_src(mgxs_lib, solver)

    # Set OpenMOC log level to reduce verbosity
    log_level = openmoc.get_log_level()
    openmoc.set_log_level('CRITICAL')

    # Initialize SPH factors
    num_groups = geometry.getNumEnergyGroups()
    num_fsrs = geometry.getNumFSRs()

    # Map FSRs to domains to compute domain-averaged fluxes
    num_domains = len(mgxs_lib.domains)
    fsrs_to_domains = np.zeros(num_fsrs)
    for fsr in range(num_fsrs):
        cell = geometry.findCellContainingFSR(fsr)

        # FIXME: Will this work with material cloning later??
        if mgxs_lib.domain_type == 'material':
            fsrs_to_domains[fsr] = cell.getFillMaterial().getId()
        else:
            fsrs_to_domains[fsr] = cell.getId()

    # Initialize array of domain-averaged fluxes and SPH factors
    openmoc_fluxes = np.zeros((num_domains, num_groups))
    sph = np.ones((num_domains, num_groups))

    # SPH Iteration
    for i in range(max_iters):

        # Run fixed source calculation
        solver.computeFlux()

        # Extract the FSR scalar fluxes and total fission rate
        fsr_fluxes = get_scalar_fluxes(solver)

        # Compute the domain-averaged flux in each energy group
        for j, domain in enumerate(mgxs_lib.domains):
            domain_fluxes = fsr_fluxes[fsrs_to_domains == domain.id, :]
            openmoc_fluxes[j, :] = np.mean(domain_fluxes, axis=0)

        # Compute SPH factors
        old_sph = np.copy(sph)
        sph = openmc_fluxes / openmoc_fluxes
        sph = np.nan_to_num(sph)
        sph[sph == 0.0] = 1.0

        # Compute SPH factor residuals
        res = np.abs((sph - old_sph) / old_sph)
        res = np.nan_to_num(res)

        # Report maximum SPH factor residual
        py_printf('RESULT', 'SPH iteration %d:\tres = %1.3e', i, res.max())

        # Update multi-group cross sections with SPH factors
        sph_mgxs_lib = _apply_sph_factors(mgxs_lib, sph, geometry)

        # Load the new MGXS library data into the OpenMOC geometry
        load_openmc_mgxs_lib(sph_mgxs_lib, geometry)

        # Check maximum SPH factor residual for convergence
        if res.max() < sph_tol:
            break

    # Warn user if SPH factors did not converge
    else:
        py_printf('WARNING', 'SPH factors did not converge')

    # Reset log level
    openmoc.set_log_level(log_level)

    return sph, sph_mgxs_lib

##
# @brief Assign fixed sources to an OpenMOC model from an OpenMC MGXS library.
# @detail this routine computes the fission production and scattering source
#         in each domain in an OpenMC MGXS library and assigns it as a fixed
#         source for an OpenMOC calculation. This is a helper routine for
#         openmoc.process.compute_sph_factors(...).
# @param mgxs_lib an openmc.mgxs.Library object
# @param solver an OpenMOC Solver object
# @return a NumPy array of the OpenMC fluxes indexed by domain and energy group
def _load_openmc_src(mgxs_lib, solver):

    import openmc

    # Retrieve dictionary of OpenMOC domains corresponding to OpenMC domains
    geometry = solver.getGeometry()
    if mgxs_lib.domain_type == 'material':
        openmoc_domains = geometry.getAllMaterials()
    else:
        openmoc_domains = geometry.getAllCells()

    # Create variables for the number of FSRs and energy groups
    num_groups = geometry.getNumEnergyGroups()
    num_domains = len(mgxs_lib.domains)
    openmc_fluxes = np.zeros((num_domains, num_groups))

    sp = openmc.StatePoint(mgxs_lib.sp_filename)
    keff = sp.k_combined[0]

    # Compute fixed sources for all domains in the MGXS library
    for i, openmc_domain in enumerate(mgxs_lib.domains):

        # Get OpenMOC domain corresponding to the OpenMOC domain
        openmoc_domain = openmoc_domains[openmc_domain.id]

        # Compute the total volume filled by this domain
        domain_volume = openmoc_domain.getVolume()
        num_instances = openmoc_domain.getNumInstances()

        # If this domain is not found in the OpenMOC geometry, ignore it
        if num_instances == 0:
            continue
        # Compute the total volume of all replicas of this Cell in the Geometry
        else:
            tot_volume = num_instances * domain_volume

        # Extract an openmc.mgxs.MGXS object for the scattering matrix
        if 'nu-scatter matrix' in mgxs_lib.mgxs_types:
            scatter = mgxs_lib.get_mgxs(openmoc_domain.getId(),
                                        'nu-scatter matrix')
        elif 'scatter matrix' in mgxs_lib.mgxs_types:
            scatter = mgxs_lib.get_mgxs(openmoc_domain.getId(),
                                        'scatter matrix')
        else:
            py_printf('ERROR', 'Unable to compute SPH factors for an OpenMC '
                               'MGXS library without scattering matrices')

        # Extract an openmc.mgxs.MGXS object for the nu-fission cross section
        if 'nu-fission' in mgxs_lib.mgxs_types:
            nu_fission = mgxs_lib.get_mgxs(openmoc_domain.getId(), 'nu-fission')
        else:
            py_printf('ERROR', 'Unable to compute SPH factors for an OpenMC '
                               'MGXS library without nu-fission cross sections')

        # Extract an openmc.mgxs.MGXS object for the chi fission spectrum
        if 'chi' in mgxs_lib.mgxs_types:
            chi = mgxs_lib.get_mgxs(openmoc_domain.getId(), 'chi')
        else:
            py_printf('ERROR', 'Unable to compute SPH factors for an OpenMC '
                               'MGXS library without nu-fission cross sections')

        # Retrieve the OpenMC volume-integrated flux for this domain from
        # the nu-fission MGXS and store it for SPH factor calculation
        flux = nu_fission.tallies['flux'].mean.flatten()
        openmc_fluxes[i, :] = np.atleast_1d(np.flipud(flux))
        openmc_fluxes[i, :] /= tot_volume

        # Extract a NumPy array for each MGXS summed across all nuclides
        scatter = scatter.get_xs(nuclides='sum')
        nu_fission = nu_fission.get_xs(nuclides='sum')
        chi = chi.get_xs(nuclides='sum')

        # Compute and store volume-averaged fission + scatter sources
        for group in range(num_groups):

            # Compute the source for this group from fission and scattering
            in_scatter = scatter[:, group] * openmc_fluxes[i, :]
            fission = (chi[group] / keff) * nu_fission * openmc_fluxes[i, :]
            source = np.sum(in_scatter) + np.sum(fission)

            # Assign the source to this domain
            if mgxs_lib.domain_type == 'material':
                solver.setFixedSourceByMaterial(openmoc_domain, group+1, source)
            else:
                solver.setFixedSourceByCell(openmoc_domain, group+1, source)

    return openmc_fluxes


##
# @brief Apply SPH factors to an OpenMC MGXS library
# @details This is a helper routine for openmoc.process.compute_sph_factors(...)
# @param mgxs_lib an openmc.mgxs.Library object
# @param sph a NumpPy array of SPH factors for each domain and energy group
# @param geometry an OpenMOC Geometry object
# @return a new openmc.mgxs.Library with SPH factors applied to all MGXS
def _apply_sph_factors(mgxs_lib, sph, geometry):

    # Retrieve dictionary of OpenMOC domains corresponding to OpenMC domains
    if mgxs_lib.domain_type == 'material':
        openmoc_domains = geometry.getAllMaterials()
    else:
        openmoc_domains = geometry.getAllCells()

    # Create a copy of the MGXS library to apply SPH factors
    sph_mgxs_lib = copy.deepcopy(mgxs_lib)

    # Loop over all domains in the MGXS library
    for i, openmc_domain in enumerate(sph_mgxs_lib.domains):
        openmoc_domain = openmoc_domains[openmc_domain.id]

        # If not a fissionable domain, then do not modify MGXS
        if not openmoc_domain.isFissionable():
            continue

        # Loop over all cross section types in the MGXS library
        for mgxs_type in sph_mgxs_lib.mgxs_types:

            # Do not update the chi fission spectrum with SPH factors
            if mgxs_type == 'chi':
                continue

            # Extract the openmc.mgxs.MGXS object from the input library
            mgxs = mgxs_lib.get_mgxs(openmoc_domain.getId(), mgxs_type)

            # Extract the openmc.mgxs.MGXS object to update with SPH factors
            sph_mgxs = sph_mgxs_lib.get_mgxs(openmoc_domain.getId(), mgxs_type)

            # Extract the OpenMC derived Tally for the MGXS
            tally = mgxs.xs_tally
            sph_tally = sph_mgxs.xs_tally
            flip_sph = np.flipud(sph[i, :])

            # If this is a scattering matrix, repeat for all outgoing groups
            if 'scatter matrix' in mgxs_type:
                flip_sph = np.repeat(flip_sph, mgxs_lib.num_groups)

            # Apply SPH factors to the MGXS in each nuclide, group
            sph_tally._mean = tally.mean * flip_sph[:, np.newaxis, np.newaxis]
            sph_tally._std_dev = tally.std_dev * flip_sph[:, np.newaxis, np.newaxis]**2

    return sph_mgxs_lib