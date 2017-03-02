import sys
import os
import copy
import collections
import hashlib

import numpy as np

import openmoc

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
    import checkvalue as cv
# For Python 3.X.X
else:
    from openmoc.log import py_printf
    import openmoc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


def _get_domain(domains, domain_spec):
    """A helper routine to find materials/cells for load_from_hdf5(...)"""

    # If domain_spec is an integer, it must be a domain ID
    if isinstance(domain_spec, int) and domain_spec in domains:
        return domains[domain_spec]

    # If domain_spec is a string, it must be a domain name
    elif isinstance(domain_spec, basestring):
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


def load_from_hdf5(filename='mgxs.h5', directory='mgxs',
                   geometry=None, domain_type='material', suffix=''):
    """This routine loads an HDF5 file of multi-group cross section data.

    The routine instantiates material with multi-group cross section data and
    returns a dictionary of each Material object keyed by its ID. An OpenMOC
    geometry may optionally be given and the routine will directly insert the
    multi-group cross sections into each material in the geometry. If a geometry
    is passed in, materials from the geometry will be used in place of those
    instantiated by this routine.

    Parameters
    ----------
    filename : str
        Filename for cross sections HDF5 file (default is 'mgxs.h5')
    directory : str
        Directory for cross sections HDF5 file (default is 'mgxs')
    geometry : openmoc.Geometry, optional
        An optional geometry populated with materials, cells, etc.
    domain_type : str
        The domain type ('material' or 'cell') upon which the cross sections
        are defined (default is 'material')
    suffix : str, optional
        An optional string suffix to index the HDF5 file beyond the assumed
        domain_type/domain_id/mgxs_type group sequence (default is '')

    Returns
    -------
    materials : dict
        A dictionary of Materials keyed by ID

    """

    cv.check_type('filename', filename, basestring)
    cv.check_type('directory', directory, basestring)
    cv.check_value('domain_type', domain_type, ('material', 'cell'))
    cv.check_type('suffix', suffix, basestring)
    if geometry:
        cv.check_type('geometry', geometry, openmoc.Geometry)

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

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    old_materials = {}
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
    for domain_spec in sorted(f[domain_type]):

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

                # If the user filled multiple Cells with the same Material,
                # the Material must be cloned for each unique Cell
                if material != None:
                    if len(domains) > geometry.getNumMaterials():
                        old_materials[material.getId()] = material
                        material = material.clone()

                # If the Cell does not contain a Material, create one for it
                else:
                    if isinstance(domain_spec, int):
                        material = openmoc.Material(id=domain_spec)
                    else:
                        # Reproducibly hash the domain name into an integer ID
                        domain_id =hashlib.md5(domain_spec.encode('utf-8'))
                        domain_id = int(domain_id.hexdigest()[:4], 16)
                        material = \
                            openmoc.Material(id=domain_id, name=domain_spec)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            if isinstance(domain_spec, int):
                material = openmoc.Material(id=domain_spec)
            else:
                # Reproducibly hash the domain name into an integer ID
                domain_id =hashlib.md5(domain_spec.encode('utf-8'))
                domain_id = int(domain_id.hexdigest()[:4], 16)
                material = openmoc.Material(id=domain_id, name=domain_spec)

        # Add material to the collection
        materials[domain_spec] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in domain_group:
            sigma = _get_numpy_array(domain_group, 'transport', suffix)
            material.setSigmaT(sigma)
            py_printf('DEBUG', 'Loaded "transport" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        elif 'total' in domain_group:
            sigma = _get_numpy_array(domain_group, 'total', suffix)
            material.setSigmaT(sigma)
            py_printf('DEBUG', 'Loaded "total" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                      '"%s %s"', domain_type, str(domain_spec))

        # Search for the fission production cross section
        if 'nu-fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-fission', suffix)
            material.setNuSigmaF(sigma)
            py_printf('DEBUG', 'Loaded "nu-fission" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                      '"%s %s"', domain_type, str(domain_spec))

        # Search for the scattering matrix cross section
        if 'nu-scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-scatter matrix', suffix)
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "nu-scatter matrix" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        elif 'scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'scatter matrix', suffix)
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "scatter matrix" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "scatter matrix" found for "%s %s"',
                      domain_type, str(domain_spec))

        # Search for chi (fission spectrum)
        if 'chi' in domain_group:
            chi = _get_numpy_array(domain_group, 'chi', suffix)
            material.setChi(chi)
            py_printf('DEBUG', 'Loaded "chi" MGXS for "%s %s"',
                      domain_type, str(domain_spec))
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %s"',
                      domain_type, str(domain_spec))

        # Search for optional cross sections
        if 'fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'fission', suffix)
            material.setSigmaF(sigma)
            py_printf('DEBUG', 'Loaded "fission" MGXS for "%s %s"',
                      domain_type, str(domain_spec))

    # Inform SWIG to garbage collect any old Materials from the Geometry
    for material_id in old_materials:
        old_materials[material_id].thisown = False

    # Return collection of materials
    return materials


def load_openmc_mgxs_lib(mgxs_lib, geometry=None):
    """This routine loads an OpenMC Library of multi-group cross section data.

    The routine instantiates materials with multi-group cross section data and
    returns a dictionary of each material keyed by its ID. An OpenMOC geometry
    may optionally be given and the routine will directly insert the multi-group
    cross sections into each material in the geometry. If a geometry is passed
    in, materials from the geometry will be used in place of those instantiated
    by this routine.

    Parameters
    ----------
    mgxs_lib : openmc.mgxs.Library
        An OpenMC multi-group cross section library library
    geometry : openmoc.Geometry, optional
        An optional geometry populated with materials, cells, etc.

    Returns
    -------
    materials : dict
        A dictionary of Materials keyed by ID

    """

    # Attempt to import openmc
    try:
        import openmc
    except ImportError:
        py_printf('ERROR', 'The OpenMC code must be installed on your system')

    cv.check_type('mgxs_lib', mgxs_lib, openmc.mgxs.Library)
    if geometry:
        cv.check_type('geometry', geometry, openmoc.Geometry)

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    old_materials = {}
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
    num_materials = 0
    for domain in mgxs_lib.domains:

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:

            # Cache the number of materials as an optimization
            if num_materials == 0:
                num_materials = geometry.getNumMaterials()

            if domain_type == 'material':
                material = _get_domain(domains, domain.id)

                # Ignore materials which cannot be found in the OpenMOC Geometry
                if material is None:
                    domain_name = domain.name.replace('%', '%%')
                    py_printf('WARNING', 'Ignoring cross sections for %s "%d" "%s"',
                              domain_type, domain.id, str(domain_name))
                    continue

            elif domain_type == 'cell':
                cell = _get_domain(domains, domain.id)

                # Ignore cells which cannot be found in the OpenMOC Geometry
                if cell is None:
                    domain_name = domain.name.replace('%', '%%')
                    py_printf('WARNING', 'Ignoring cross sections for %s "%d" "%s"',
                              domain_type, domain.id, str(domain_name))
                    continue
                else:
                    material = cell.getFillMaterial()

                # If the user filled multiple Cells with the same Material,
                # the Material must be cloned for each unique Cell
                if material != None:
                    if len(domains) > num_materials:
                        old_materials[material.getId()] = material
                        material = material.clone()

                # If the Cell does not contain a Material, create one for it
                else:
                    material = openmoc.Material(id=domain.id)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            material = openmoc.Material(id=domain.id)

        domain_name = domain.name.replace('%', '%%')
        py_printf('INFO', 'Importing cross sections for %s "%d" "%s"',
                  domain_type, domain.id, str(domain_name))

        # Add material to the collection
        materials[domain.id] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'transport')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
            py_printf('DEBUG', 'Loaded "transport" MGXS for "%s %d"',
                      domain_type, domain.id)
        elif 'nu-transport' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-transport')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
            py_printf('DEBUG', 'Loaded "nu-transport" MGXS for "%s %d"',
                      domain_type, domain.id)
        elif 'total' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'total')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
            py_printf('DEBUG', 'Loaded "total" MGXS for "%s %d"',
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                      '"%s %d"', domain_type, domain.id)

        # Search for the fission production cross section
        if 'nu-fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setNuSigmaF(sigma)
            py_printf('DEBUG', 'Loaded "nu-fission" MGXS for "%s %d"',
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                      '"%s %d"', domain_type, domain.id)

        # Search for the scattering matrix cross section
        if 'consistent nu-scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'consistent nu-scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "consistent nu-scatter matrix" MGXS for "%s %d"',
                      domain_type, domain.id)
        elif 'nu-scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "nu-scatter matrix" MGXS for "%s %d"',
                      domain_type, domain.id)
        elif 'consistent scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'consistent scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "consistent scatter matrix" MGXS for "%s %d"',
                      domain_type, domain.id)
        elif 'scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
            py_printf('DEBUG', 'Loaded "scatter matrix" MGXS for "%s %d"',
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "scatter matrix" or "nu-scatter matrix" '
                      'found for "%s %d"', domain_type, domain.id)

        # Search for chi (fission spectrum)
        if 'chi' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'chi')
            chi = mgxs.get_xs(nuclides='sum')
            material.setChi(chi)
            py_printf('DEBUG', 'Loaded "chi" MGXS for "%s %d"',
                      domain_type, domain.id)
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %d"',
                      domain_type, domain.id)

        # Search for optional cross sections
        if 'fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaF(sigma)
            py_printf('DEBUG', 'Loaded "fission" MGXS for "%s %d"',
                      domain_type, domain.id)

    # Inform SWIG to garbage collect any old Materials from the Geometry
    for material_id in old_materials:
        old_materials[material_id].thisown = False

    # Return collection of materials
    return materials


def compute_sph_factors(mgxs_lib, max_sph_iters=30, sph_tol=1E-5,
                        fix_src_tol=1E-5, num_azim=4, azim_spacing=0.1,
                        zcoord=0.0, num_threads=1, throttle_output=True):
    """Compute SPH factors for an OpenMC multi-group cross section library.

    This routine coputes SuPerHomogenisation (SPH) factors for an OpenMC MGXS
    library. The SPH scheme is outlined by Alain Hebert in the following paper:

        Hebert, A., "A Consistent Technique for the Pin-by-Pin
        Homogenization of a Pressurized Water Reactor Assembly."
        Nuclear Science and Engineering, 113 (3), pp. 227-233, 1993.

    The SPH factors are needed to preserve reaction rates in heterogeneous
    geometries. The energy condensation process leads to a bias between
    ultrafine and coarse energy group calculations. This bias is a result of the
    use of scalar flux-weighting to compute MGXS without properly accounting for
    angular-dependence of the flux.

    Parameters
    ----------
    mgxs_lib : openmc.mgxs.Library
        An OpenMC multi-group cross section library
    max_sph_iters : Integral
        The maximum number of SPH iterations (default is 30)
    sph_tol : Real
        The tolerance on the SPH factor convergence (default is 1E-5)
    fix_src_tol : Real
        The tolerance on the MOC fixed source calculations (default is 1E-5_
    num_azim : Integral
        The number of azimuthal angles (default is 4)
    azim_spacing : Real
        The track spacing (default is 0.1 centimeters)
    zcoord : Real
        The coordinate on the z-axis (default is 0.)
    num_threads : Real
        The number of OpenMP threads (default is 1)
    throttle_output : bool
        Whether to suppress output from fixed source calculations (default is True)

    Returns
    -------
    fsrs_to_sph : numpy.ndarray of Real
        A NumPy array of SPH factors indexed by FSR and energy group
    sph_mgxs_lib : openmc.mgxs.Library
        A OpenMC MGXS library with the SPH factors applied to each MGXS
    sph_to_fsrs_indices : numpy.ndarray of Integral
        A NumPy array of all fissionable FSRs to which SPH factors were applied

    """

    import openmc.mgxs

    cv.check_type('mgxs_lib', mgxs_lib, openmc.mgxs.Library)

    from openmc.openmoc_compatible import get_openmoc_geometry

    # For Python 2.X.X
    if sys.version_info[0] == 2:
        from process import get_scalar_fluxes
    # For Python 3.X.X
    else:
        from openmoc.process import get_scalar_fluxes

    py_printf('NORMAL', 'Computing SPH factors...')

    # Create an OpenMOC Geometry from the OpenMC Geometry
    geometry = get_openmoc_geometry(mgxs_lib.geometry)

    # Load the MGXS library data into the OpenMOC geometry
    load_openmc_mgxs_lib(mgxs_lib, geometry)

    # Initialize an OpenMOC TrackGenerator
    track_generator = openmoc.TrackGenerator(geometry, num_azim, azim_spacing)
    track_generator.setZCoord(zcoord)
    track_generator.generateTracks()

    # Initialize an OpenMOC Solver
    solver = openmoc.CPUSolver(track_generator)
    solver.setConvergenceThreshold(fix_src_tol)
    solver.setNumThreads(num_threads)
    openmc_fluxes = _load_openmc_src(mgxs_lib, solver)

    # Initialize SPH factors
    num_groups = geometry.getNumEnergyGroups()
    num_fsrs = geometry.getNumFSRs()

    # Map FSRs to domains (and vice versa) to compute domain-averaged fluxes
    fsrs_to_domains = np.zeros(num_fsrs)
    domains_to_fsrs = collections.defaultdict(list)
    sph_to_fsr_indices = []

    for fsr in range(num_fsrs):
        cell = geometry.findCellContainingFSR(fsr)

        if mgxs_lib.domain_type == 'material':
            domain = cell.getFillMaterial()
        else:
            domain = cell

        fsrs_to_domains[fsr] = domain.getId()
        domains_to_fsrs[domain.getId()].append(fsr)

        if domain.isFissionable():
            sph_to_fsr_indices.append(fsr)

    # Get all OpenMOC domains
    if mgxs_lib.domain_type == 'material':
        openmoc_domains = geometry.getAllMaterials()
    elif mgxs_lib.domain_type == 'cell':
        openmoc_domains = geometry.getAllMaterialCells()
    else:
        py_printf('ERROR', 'SPH factors cannot be applied for an OpenMC MGXS '
                  'library of domain type %s', mgxs_lib.domain_type)

    # Build a list of indices into the SPH array for fissionable domains
    sph_to_domain_indices = []
    for i, openmc_domain in enumerate(mgxs_lib.domains):
        if openmc_domain.id in openmoc_domains:
            openmoc_domain = openmoc_domains[openmc_domain.id]
            if openmoc_domain.isFissionable():
                sph_to_domain_indices.append(i)

    py_printf('NORMAL', 'Computing SPH factors for %d "%s" domains',
               len(sph_to_domain_indices), mgxs_lib.domain_type)

    # Initialize array of domain-averaged fluxes and SPH factors
    num_domains = len(mgxs_lib.domains)
    openmoc_fluxes = np.zeros((num_domains, num_groups))
    sph = np.ones((num_domains, num_groups))

    # Store starting verbosity log level
    log_level = openmoc.get_log_level()

    # SPH iteration loop
    for i in range(max_sph_iters):

        # Run fixed source calculation with suppressed output
        if throttle_output:
            openmoc.set_log_level('WARNING')

        # Fixed source calculation
        solver.computeFlux()

        # Restore log output level
        if throttle_output:
            openmoc.set_log_level(log_level)

        # Extract the FSR scalar fluxes
        fsr_fluxes = get_scalar_fluxes(solver)

        # Compute the domain-averaged flux in each energy group
        for j, openmc_domain in enumerate(mgxs_lib.domains):
            domain_fluxes = fsr_fluxes[fsrs_to_domains == openmc_domain.id, :]
            openmoc_fluxes[j, :] = np.mean(domain_fluxes, axis=0)

        # Compute SPH factors
        old_sph = np.copy(sph)
        sph = openmc_fluxes / openmoc_fluxes
        sph = np.nan_to_num(sph)
        sph[sph == 0.0] = 1.0

        # Compute SPH factor residuals
        res = np.abs((sph - old_sph) / old_sph)
        res = np.nan_to_num(res)

        # Extract residuals for fissionable domains only
        res = res[sph_to_domain_indices, :]

        # Report maximum SPH factor residual
        py_printf('NORMAL', 'SPH Iteration %d:\tres = %1.3e', i, res.max())

        # Create a new MGXS library with cross sections updated by SPH factors
        sph_mgxs_lib = _apply_sph_factors(mgxs_lib, geometry, sph)

        # Load the new MGXS library data into the OpenMOC geometry
        load_openmc_mgxs_lib(sph_mgxs_lib, geometry)

        # Check max SPH factor residual for this domain for convergence
        if res.max() < sph_tol:
            break

    # Warn user if SPH factors did not converge
    else:
        py_printf('WARNING', 'SPH factors did not converge')

    # Collect SPH factors for each FSR, energy group
    fsrs_to_sph = np.ones((num_fsrs, num_groups), dtype=np.float)
    for i, openmc_domain in enumerate(mgxs_lib.domains):
        if openmc_domain.id in openmoc_domains:
            openmoc_domain = openmoc_domains[openmc_domain.id]
            if openmoc_domain.isFissionable():
                fsr_ids = domains_to_fsrs[openmc_domain.id]
                fsrs_to_sph[fsr_ids,:] = sph[i,:]

    return fsrs_to_sph, sph_mgxs_lib, np.array(sph_to_fsr_indices)


def _load_openmc_src(mgxs_lib, solver):
    """Assign fixed sources to an OpenMOC model from an OpenMC MGXS library.

    This routine computes the fission production and scattering source in
    each domain in an OpenMC MGXS library and assigns it as a fixed source
    for an OpenMOC calculation. This is a helper routine for the
    compute_sph_factors(...) routine.

    Parameters
    ----------
    mgxs_lib : openmc.mgxs.Library object
        An OpenMC multi-group cross section library
    solver : openmoc.Solver
        An OpenMOC solver into which to load the fixed sources

    Returns
    -------
    openmc_fluxes : numpy.ndarray of Real
        A NumPy array of the OpenMC fluxes indexed by domain and energy group

    """

    # Retrieve dictionary of OpenMOC domains corresponding to OpenMC domains
    geometry = solver.getGeometry()
    if mgxs_lib.domain_type == 'material':
        openmoc_domains = geometry.getAllMaterials()
    else:
        openmoc_domains = geometry.getAllCells()

    # Create variables for the number of domains and energy groups
    num_groups = geometry.getNumEnergyGroups()
    num_domains = len(mgxs_lib.domains)
    openmc_fluxes = np.zeros((num_domains, num_groups))
    keff = mgxs_lib.keff

    # Create mapping of FSRs-to-domains to optimize fixed source setup
    domains_to_fsrs = collections.defaultdict(list)
    for fsr_id in range(geometry.getNumFSRs()):
        cell = geometry.findCellContainingFSR(fsr_id)
        if mgxs_lib.domain_type == 'material':
            domain = cell.getFillMaterial()
        else:
            domain = cell
        domains_to_fsrs[domain.getId()].append(fsr_id)

    # Compute fixed sources for all domains in the MGXS library
    for i, openmc_domain in enumerate(mgxs_lib.domains):

        # Ignore domains which cannot be found in the OpenMOC Geometry
        if openmc_domain.id not in openmoc_domains:
            continue

        # Get OpenMOC domain corresponding to the OpenMC domain
        openmoc_domain = openmoc_domains[openmc_domain.id]

        # If this domain is not found in the OpenMOC geometry, ignore it
        if openmoc_domain.getNumInstances() == 0:
            continue

        # Compute the total volume filled by this domain throughout the geometry
        tot_volume = openmoc_domain.getVolume()

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
                      'MGXS library without chi fission spectrum')

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
            fsr_ids = domains_to_fsrs[openmoc_domain.getId()]
            for fsr_id in fsr_ids:
                solver.setFixedSourceByFSR(fsr_id, group+1, source)

    return openmc_fluxes


def _apply_sph_factors(mgxs_lib, geometry, sph):
    """Apply SPH factors to an OpenMC MGXS library.

    This is a helper routine for the compute_sph_factors(...) routine.

    Parameters
    ----------
    mgxs_lib : openmc.mgxs.Library
        An OpenMC multi-group cross section library
    geometry : openmoc.Geometry
        An OpenMOC geometry
    sph : numpy.ndarray of Real
        A NumpPy array of SPH factors for each domain and energy group

    Returns
    -------
    sph_mgxs_lib : openmc.mgxs.Library
        A new OpenMC MGXS library with SPH factors applied to all MGXS

    """

    # Create a copy of the MGXS library to apply SPH factors
    sph_mgxs_lib = copy.deepcopy(mgxs_lib)

    # Get all OpenMOC domains
    if mgxs_lib.domain_type == 'material':
        openmoc_domains = geometry.getAllMaterials()
    else:
        openmoc_domains = geometry.getAllMaterialCells()

    # Loop over all domains
    for i, openmc_domain in enumerate(mgxs_lib.domains):

        # Ignore domains which cannot be found in the OpenMOC Geometry
        if openmc_domain.id not in openmoc_domains:
            continue

        # Get OpenMOC domain corresponding to the OpenMC domain
        openmoc_domain = openmoc_domains[openmc_domain.id]

        # Ignore non-fissionable domains
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
            flip_sph = np.flipud(sph[i,:])

            # If this is a scattering matrix, repeat for all outgoing groups
            if 'scatter matrix' in mgxs_type:
                flip_sph = np.repeat(flip_sph, mgxs_lib.num_groups)

            # Apply SPH factors to the MGXS in each nuclide, group
            sph_tally._mean = tally.mean * flip_sph[:, np.newaxis, np.newaxis]
            sph_tally._std_dev = \
                tally.std_dev * flip_sph[:, np.newaxis, np.newaxis]**2

    return sph_mgxs_lib
