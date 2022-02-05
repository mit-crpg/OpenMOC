import os
import re
import mmap
import sys
import math
import datetime
import operator
from numbers import Integral, Real
import pickle

import numpy as np
import h5py

import openmoc

# For Python 2.X.X
if (sys.version_info[0] == 2):
    from log import *
    import checkvalue as cv
# For Python 3.X.X
else:
    from openmoc.log import *
    import openmoc.checkvalue as cv

if sys.version_info >= (3, 3):
    from collections.abc import Iterable
else:
    from collections import Iterable

# Store viable OpenMOC solver types for type checking
solver_types = (openmoc.Solver,)
# Store implemented reaction types for value checking
rxn_types = ('flux', 'total', 'scatter', 'fission', 'nu-fission')
try:
    # Try to import OpenMOC's CUDA module
    if (sys.version_info[0] == 2):
        from cuda import GPUSolver
    else:
        from openmoc.cuda import GPUSolver
    solver_types += (GPUSolver,)
except ImportError:
    pass

if sys.version_info[0] >= 3:
    basestring = str


def get_scalar_fluxes(solver, fsrs='all', groups='all'):
    """Return an array of scalar fluxes in one or more FSRs and groups.

    This routine builds a 2D NumPy array indexed by FSR and energy group for
    the corresponding scalar fluxes. The fluxes are organized in the array in
    order of increasing FSR and enery group if 'all' FSRs or energy groups are
    requested (the default). If the user requests fluxes for specific FSRs or
    energy groups, then the fluxes are returned in the order in which the FSRs
    and groups are enumerated in the associated paramters.

    Parameters
    ----------
    solver : openmoc.Solver
        The solver used to compute the flux
    fsrs : Iterable of Integral or 'all'
        A collection of integer FSR IDs or 'all' (default)
    groups : Iterable of Integral or 'all'
        A collection of integer energy groups or 'all' (default)

    Returns
    -------
    fluxes : ndarray
        The scalar fluxes indexed by FSR ID and energy group. Note that the
        energy group index starts at 0 rather than 1 for the highest energy
        in accordance with Python's 0-based indexing.

    """

    global solver_types
    cv.check_type('solver', solver, solver_types)

    if isinstance('fsrs', basestring):
        cv.check_value('fsrs', fsrs, 'all')
    else:
        cv.check_type('fsrs', Iterable, Integral)

    if isinstance('groups', basestring):
        cv.check_value('groups', fsrs, 'all')
    else:
        cv.check_type('groups', Iterable, Integral)

    # Extract all of the FSR scalar fluxes
    if groups == 'all' and fsrs == 'all':
        num_fsrs = int(solver.getGeometry().getNumTotalFSRs())
        num_groups = solver.getGeometry().getNumEnergyGroups()
        num_fluxes = num_groups * num_fsrs
        fluxes = solver.getFluxes(num_fluxes)
        fluxes = np.reshape(fluxes, (num_fsrs, num_groups))
        return fluxes

    # Build a list of FSRs to iterate over
    if fsrs == 'all':
        num_fsrs = int(solver.getGeometry().getNumTotalFSRs())
        fsrs = np.arange(num_fsrs)
    else:
        num_fsrs = len(fsrs)

    # Build a list of enery groups to iterate over
    if groups == 'all':
        num_groups = solver.getGeometry().getNumEnergyGroups()
        groups = np.arange(num_groups) + 1
    else:
        num_groups = len(groups)

    # Extract some of the FSR scalar fluxes
    fluxes = np.zeros((num_fsrs, num_groups))
    for fsr in fsrs:
        for group in groups:
            fluxes[fsr, group-1] = solver.getFlux(int(fsr), int(group))

    return fluxes


def compute_fission_rates(solver, use_hdf5=False):
    """Computes the fission rate in each FSR.

    This method combines the rates based on their hierarchical universe/lattice
    structure. The fission rates are then exported to a binary HDF5 or Python
    pickle file.

    This routine is intended to be called by the user in Python to compute
    fission rates. Typically, the fission rates will represent pin powers. The
    routine either exports fission rates to an HDF5 binary file or pickle file
    with each fission rate being indexed by a string representing the
    universe/lattice hierarchy.

    Parameters
    ----------
    solver : openmoc.Solver
        The solver used to compute the flux
    use_hdf5 : bool
        Whether or not to export fission rates to an HDF5 file

    Examples
    --------
    This routine may be called from a Python script as follows:

        >>> compute_fission_rates(solver, use_hdf5=True)

    """

    global solver_types
    cv.check_type('solver', solver, solver_types)
    cv.check_type('use_hdf5', use_hdf5, bool)

    # Make directory if it does not exist
    directory = openmoc.get_output_directory() + '/fission-rates/'
    filename = 'fission-rates'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Get geometry
    geometry = solver.getGeometry()

    # Compute the volume-weighted fission rates for each FSR
    fsr_fission_rates = \
        solver.computeFSRFissionRates(int(geometry.getNumTotalFSRs()))

    # Initialize fission rates dictionary
    fission_rates_sum = {}

    # Loop over FSRs and populate fission rates dictionary
    for fsr in range(int(geometry.getNumTotalFSRs())):

        if geometry.findFSRMaterial(fsr).isFissionable():

            # Get the linked list of LocalCoords
            point = geometry.getFSRPoint(fsr)
            coords = openmoc.LocalCoords(point.getX(), point.getY(), point.getZ())
            coords.setUniverse(geometry.getRootUniverse())
            geometry.findCellContainingCoords(coords)
            coords = coords.getHighestLevel().getNext()

            # initialize dictionary key
            key = 'UNIV = 0 : '

            # Parse through the linked list and create fsr key.
            # If lowest level sub dictionary already exists, then increment
            # fission rate; otherwise, set the fission rate.
            while True:
                if coords.getType() == openmoc.LAT:
                    key += 'LAT = ' + str(coords.getLattice().getId()) + ' (' + \
                           str(coords.getLatticeX()) + ', ' + \
                           str(coords.getLatticeY()) + ', ' + \
                           str(coords.getLatticeZ()) + ') : '
                else:
                    key += 'UNIV = ' + str(coords.getUniverse().getId()) + ' : '

                # Remove trailing ' : ' on end of key if at last univ/lat
                if coords.getNext() is None:
                    key = key[:-3]
                    break
                else:
                    coords = coords.getNext()

            # Increment or set fission rate
            if key in fission_rates_sum:
                fission_rates_sum[key] += fsr_fission_rates[fsr]
            else:
                fission_rates_sum[key] = fsr_fission_rates[fsr]

    # Write the fission rates to the HDF5 file
    if use_hdf5:
        f = h5py.File(directory + filename + '.h5', 'w')
        fission_rates_group = f.create_group('fission-rates')
        for key, value in fission_rates_sum.items():
            fission_rates_group.attrs[key] = value
        f.close()

    # Pickle the fission rates to a file
    else:
        pickle.dump(fission_rates_sum, open(directory + filename + '.pkl', 'wb'))


def get_sigma_by_group(material, rxn_type, g):
    """Return the cross section for a given reaction

    Parameters
    ----------
    material : openmoc.Material
        Material that we are getting the cross sections from
    rxn_type : {'flux', 'total', 'scatter',
                'fission', 'nu-fission'}
        Reaction type we are loading the cross section for
    g : integer
        Energy group index (starts at 0).

    Returns
    -------
    sigma : float
        The cross section for energy group `g' of reaction `rxn_type'.
        For scatter, this includes self-scatter and outscatter.

    """
    global rxn_types
    cv.check_value('rxn_type', rxn_type, rxn_types)

    # Energy groups start at 1 in OpenMOC
    if rxn_type == "total":
        return material.getSigmaTByGroup(g + 1)
    elif rxn_type == "fission":
        return material.getSigmaFByGroup(g + 1)
    elif rxn_type == "nu-fission":
        return material.getNuSigmaFByGroup(g + 1)
    elif rxn_type == "scatter":
        scatter = 0.
        for gprime in range(material.getNumEnergyGroups()):
            scatter += material.getSigmaSByGroup(g + 1, gprime + 1)
        return scatter
    else:
        # Flux
        return 1.


def store_simulation_state(solver, fluxes=False, sources=False,
                           fission_rates=False, use_hdf5=False,
                           filename='simulation-state',
                           directory='simulation-states',
                           append=True, note=''):
    """Store all of the data for an OpenMOC simulation to a binary file for
    downstream data processing.

    This routine may be used to store the following:

        * type of Solver used
        * floating point precision
        * exponential evaluation method
        * number of FSRs
        * number of materials
        * number of energy groups
        * number of azimuthal angles
        * number of polar angles
        * track spacing
        * number of tracks
        * number of track segments
        * number of source iterations
        * source convergence tolerance
        * converged $k_{eff}$
        * total runtime [seconds]
        * number of OpenMP or CUDA threads

    In addition, the routine can optionally store the FSR scalar fluxes, FSR
    sources, and pin and assembly fission rates.

    The routine may export the simulation data to either an HDF5 or a Python
    pickle binary file. Users may tell the routine to either create a new binary
    output file, or append to an existing file using a timestamp to record
    multiple simulation states to the same file.

    Parameters
    ----------
    solver : openmoc.Solver
        The solver used to compute the flux
    fluxes : bool
        Whether to store FSR scalar fluxes (False by default)
    sources : bool
        Whether to store FSR sources (False by default)
    fission_rates : bool
        Whether to store fission rates (False by default)
    use_hdf5 : bool
        Whether to export to HDF5 (True by default) or Python pickle file
    filename : str
        The filename to use (default is 'simulation-state.h5')
    directory : str
        The directory to use (default is 'simulation-states')
    append : bool
        Append to existing file or create new one (False by default)
    note : str, optional
        An optional string note to include in state file

    Examples
    --------
    This routine may be called from Python as follows:

        >>> store_simulation_state(solver, fluxes=True, source=True, \
                                   fission_rates=True, use_hdf5=True)

    See Also
    --------
    restore_simulation_state(...)

    """

    global solver_types
    cv.check_type('solver', solver, solver_types)
    cv.check_type('fluxes', fluxes, bool)
    cv.check_type('sources', sources, bool)
    cv.check_type('fission_rates', fission_rates, bool)
    cv.check_type('use_hdf5', use_hdf5, bool)
    cv.check_type('filename', filename, basestring)
    cv.check_type('directory', directory, basestring)
    cv.check_type('append', append, bool)
    cv.check_type('note', note, basestring)

    # Make directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Get the day and time to construct the appropriate groups in the file
    time = datetime.datetime.now()
    year = time.year
    month = time.month
    day = time.day
    hr = time.hour
    mins = time.minute
    sec = time.second

    # Determine the Solver type
    solver_type = ''

    if 'CPUSolver' in str(solver.__class__):
        solver_type = 'CPUSolver'
    elif 'VectorizedSolver' in str(solver.__class__):
        solver_type = 'VectorizedSolver'
    elif 'GPUSolver' in str(solver.__class__):
        solver_type = 'GPUSolver'

    # Determine the floating point precision level
    if solver.isUsingDoublePrecision():
        precision = 'double'
    else:
        precision = 'single'

    # Determine whether we are using the exponential
    # linear interpolation for exponential evaluations
    if solver.isUsingExponentialInterpolation():
        method = 'linear interpolation'
    else:
        method = 'exp intrinsic'

    # Determine whether the Solver has initialized Coarse Mesh Finite
    # Difference Acceleration (CMFD)
    if solver.getGeometry().getCmfd() is not None:
        cmfd = True
    else:
        cmfd = False

    # Get the Geometry and TrackGenerator from the solver
    geometry = solver.getGeometry()
    track_generator = solver.getTrackGenerator()

    # Retrieve useful data from the Solver, Geometry and TrackGenerator
    num_FSRs = int(geometry.getNumTotalFSRs())
    num_materials = geometry.getNumMaterials()
    num_groups = geometry.getNumEnergyGroups()
    zcoord = track_generator.getZCoord()
    num_tracks = track_generator.getNumTracks()
    num_segments = track_generator.getNumSegments()
    spacing = track_generator.getDesiredAzimSpacing()
    num_azim = track_generator.getNumAzim()
    num_polar = solver.getNumPolarAngles()
    num_iters = solver.getNumIterations()
    thresh = solver.getConvergenceThreshold()
    tot_time = solver.getTotalTime()
    keff = solver.getKeff()

    if solver_type == 'GPUSolver':
        num_threads = solver.getNumThreadsPerBlock()
        num_blocks = solver.getNumThreadBlocks()
    else:
        num_threads = solver.getNumThreads()

    # If the user requested to store the FSR fluxes
    if fluxes:

        # Allocate array
        scalar_fluxes = np.zeros((num_FSRs, num_groups))

        # Get the scalar flux for each FSR and energy group
        for i in range(num_FSRs):
            for j in range(num_groups):
                scalar_fluxes[i,j] = solver.getFlux(i,j+1)

    # If the user requested to store the FSR sources
    if sources:

        # Allocate array
        sources_array = np.zeros((num_FSRs, num_groups))

        # Get the scalar flux for each FSR and energy group
        for i in range(num_FSRs):
            for j in range(num_groups):
                sources_array[i,j] = solver.getFSRSource(i,j+1)

    # If using HDF5
    if use_hdf5:
        if append:
            f = h5py.File(directory + '/' + filename + '.h5', 'a')
        else:
            f = h5py.File(directory + '/' + filename + '.h5', 'w')

        # Create groups for the day in the HDF5 file
        day_key = '{0:02}-{1:02}-{2:02}'.format(month, day, year)
        day_group = f.require_group(day_key)

        # Create group for the time - use counter in case two simulations
        # write simulation state at the exact same hour,minute, and second
        time_key = '{0:02}:{1:02}:{2:02}'.format(hr, mins, sec)
        counter = 0
        while time_key in day_group.keys():
            time_key = '{0:02}:{1:02}:{2:02}-{3}'.format(hr, mins, sec, counter)
            counter += 1

        time_group = day_group.require_group(time_key)

        # Store a note for this simulation state
        if not note == '':
            time_group.attrs['note'] = note

        # Store simulation data to the HDF5 file
        time_group.create_dataset('solver type', data=solver_type)
        time_group.create_dataset('# FSRs', data=num_FSRs)
        time_group.create_dataset('# materials', data=num_materials)
        time_group.create_dataset('# energy groups', data=num_groups)
        time_group.create_dataset('z coord', data=zcoord)
        time_group.create_dataset('# tracks', data=num_tracks)
        time_group.create_dataset('# segments', data=num_segments)
        time_group.create_dataset('track spacing [cm]', data=spacing)
        time_group.create_dataset('# azimuthal angles', data=num_azim)
        time_group.create_dataset('# polar angles', data=num_polar)
        time_group.create_dataset('# iterations', data=num_iters)
        time_group.create_dataset('convergence threshold', data=thresh)
        time_group.create_dataset('exponential', data=method)
        time_group.create_dataset('floating point', data=precision)
        time_group.create_dataset('CMFD', data=cmfd)
        time_group.create_dataset('time [sec]', data=tot_time)
        time_group.create_dataset('keff', data=keff)

        if solver_type == 'GPUSolver':
            time_group.create_dataset('# threads per block', data=num_threads)
            time_group.create_dataset('# thread blocks', data=num_blocks)
        else:
            time_group.create_dataset('# threads', data=num_threads)

        if fluxes:
            time_group.create_dataset('FSR scalar fluxes', data=scalar_fluxes)

        if sources:
            time_group.create_dataset('FSR sources', data=sources_array)

        if fission_rates:
            compute_fission_rates(solver, use_hdf5=True)
            fission_rates_file = h5py.File('fission-rates/fission-rates.h5', 'r')
            f.copy(fission_rates_file, time_group, name='fission-rates')
            fission_rates_file.close()

        # Close the HDF5 file
        f.close()

    # If not using HDF5, we are pickling all of the data
    else:
        filename = directory + '/' + filename + '.pkl'
        if os.path.exists(filename) and append:
            sim_states = pickle.load(open(filename, 'rb'))
        else:
            sim_states = {}

        # Create strings for the day and time
        day = str(month).zfill(2)+'-'+str(day).zfill(2)+'-'+str(year)
        time = str(hr).zfill(2)+':'+str(mins).zfill(2)+':'+str(sec).zfill(2)

        # Create dictionaries for this day and time within the pickled file
        if not day in sim_states.keys():
            sim_states[day] = {}

        sim_states[day][time] = {}
        state = sim_states[day][time]

        # Store a note for this simulation state
        if not note == '':
            state['note'] = note

        # Store simulation data to a Python dictionary
        state['solver type'] = solver_type
        state['# FSRs'] = num_FSRs
        state['# materials'] = num_materials
        state['# energy groups'] = num_groups
        state['z coord'] = zcoord
        state['# tracks'] = num_tracks
        state['# segments'] = num_segments
        state['track spacing [cm]'] = spacing
        state['# azimuthal angles'] = num_azim
        state['# polar angles'] = num_polar
        state['# iterations'] = num_iters
        state['convergence threshold'] = thresh
        state['exponential'] = method
        state['floating point'] = precision
        state['CMFD'] = cmfd
        state['time [sec]'] = tot_time
        state['keff'] = keff

        if solver_type == 'GPUSolver':
            state['# threads per block'] = num_threads
            state['# thread blocks'] = num_blocks
        else:
            state['# threads'] = num_threads

        if fluxes:
            state['FSR scalar fluxes'] = scalar_fluxes

        if sources:
            state['FSR sources'] = sources_array

        if fission_rates:
            compute_fission_rates(solver, False)
            state['fission-rates'] = \
                pickle.load(open('fission-rates/fission-rates.pkl', 'rb'))

        # Pickle the simulation states to a file
        pickle.dump(sim_states, open(filename, 'wb'))

        # Pickle the simulation states to a file
        pickle.dump(sim_states, open(filename, 'wb'))


def restore_simulation_state(filename='simulation-state.h5',
                             directory='simulation-states'):
    """Restore all of the data for an OpenMOC simulation from a binary file for
    downstream data processing to a Python dictionary.

    This routine may import the simulation state from either an HDF5 or a Python
    pickle binary file created by the store_simulation_state(...) method. The
    method may be used to restore the following information:

        * type of Solver used
        * floating point precision
        * exponential evaluation method
        * number of FSRs
        * number of materials
        * number of energy groups
        * number of azimuthal angles
        * number of polar angles
        * track spacing
        * number of tracks
        * number of track segments
        * number of source iterations
        * source convergence tolerance
        * converged $k_{eff}$
        * total runtime [seconds]
        * number of OpenMP or CUDA threads

          Note: If the fission rates were stored in a hdf5 binary file,
          they are not restored and returned in this method.

    Paramters
    ---------
    filename : str
        The simulation state filename string
    directory : str
        The directory where to find the simulation state file

    Returns
    -------
    states : dict
        The dictionary of key/value pairs for simulation state data

    Examples
    --------
    This method may be called from Python as follows:

        >>> restore_simulation_state(filename='simulation-state-v1.3.h5')

    See Also
    --------
    store_simulation_state(...)

    """

    cv.check_type('filename', filename, basestring)
    cv.check_type('directory', directory, basestring)

    filename = directory + '/' + filename
    if not os.path.isfile(filename):
        py_printf('ERROR', 'Unable restore simulation state since "{0}" ' + \
                  'is not an existing simulation state file'.format(filename))

    # If using HDF5
    if '.h5' in filename or '.hdf5' in filename:

        import h5py

        # Create a file handle
        f = h5py.File(filename, 'r')

        states = {}

        # Loop over all simulation state timestamps by day
        for day in f.keys():

            # Create sub-dictionary for this day
            states[day] = {}

            # Loop over all simulation state timestamps by time of day
            for time in f[day]:

                # Create sub-dictionary for this simulation state
                dataset = f[day][time]
                states[day][time] = {}
                state = states[day][time]

                # Extract simulation state data
                solver_type = str(dataset['solver type'])
                num_FSRs = int(dataset['# FSRs'][...])
                num_materials = int(dataset['# materials'][...])
                num_tracks = int(dataset['# tracks'][...])
                num_segments = int(dataset['# segments'][...])
                spacing = int(dataset['track spacing [cm]'][...])
                num_azim = int(dataset['# azimuthal angles'][...])
                num_polar = int(dataset['# polar angles'][...])
                num_iters = int(dataset['# iterations'][...])
                thresh = float(dataset['convergence threshold'][...])
                method = str(dataset['exponential'][...])
                precision = str(dataset['floating point'][...])
                cmfd = str(dataset['CMFD'][...])
                time = float(dataset['time [sec]'][...])
                keff =  float(dataset['keff'][...])

                # Store simulation state data in sub-dictionary
                state['solver type'] = solver_type
                state['# FSRs'] = num_FSRs
                state['# materials'] = num_materials
                state['# tracks'] = num_tracks
                state['# segments'] = num_segments
                state['track spacing [cm]'] = spacing
                state['# azimuthal angles'] = num_azim
                state['# polar angles'] = num_polar
                state['# iterations'] = num_iters
                state['convergence threshold'] = thresh
                state['exponential'] = method
                state['floating point'] = precision
                state['CMFD'] = cmfd
                state['time [sec]'] = time
                state['keff'] = keff

                if solver_type == 'GPUSolver':
                    state['# threads per block'] = \
                        int(dataset['# threads per block'])
                    state['# thread blocks'] = int(dataset['# thread blocks'])
                else:
                    state['# threads'] = int(dataset['# threads'][...])
                if 'FSR scalar fluxes' in dataset:
                    state['FSR scalar fluxes'] = \
                        dataset['FSR scalar fluxes'][...]
                if 'FSR sources' in dataset:
                    state['FSR sources'] = dataset['FSR sources'][...]
                if 'note' in dataset:
                    state['note'] = str(dataset['note'])
                if 'fission-rates' in dataset:
                    py_printf('WARNING', 'The restore_simulation_state(...)' +
                              'method does not yet support fission rates')

        return states

    # If using a Python pickled file
    elif '.pkl' in filename:
        states = pickle.load(open(filename, 'rb'))
        return states

    # If file does not have a recognizable extension
    else:
        py_printf('WARNING', 'Unable to restore the simulation states file %s' +
                  ' since it does not have a supported file extension. Only ' +
                  '*.h5, *.hdf5, and *.pkl files are supported', filename)
        return {}


def parse_convergence_data(filename, directory=''):
    """Parse an OpenMOC log file to obtain a simulation's convergence data.

    This method compiles the eigenvalue and source residuals from each iteration
    of an OpenMOC simulation. This data is inserted into a Python dictionary
    under the key names 'eigenvalues' and 'residuals', along with an integer
    '# iters', and returned to the user.

    Parameters
    ----------
    filename : str
        The OpenMOC log filename string
    directory : str
        The directory where to find the log file

    Returns
    -------
    convergence_data : dict
        A Python dictionary of key/value pairs for convergence data

    Examples
    --------
    This method may be called from Python as follows:

        >>> parse_convergence_data(filename='openmoc-XX-XX-XXXX--XX:XX:XX.log')

    """

    cv.check_type('filename', filename, basestring)
    cv.check_type('directory', directory, basestring)

    # If the user specified a directory
    if len(directory) > 0:
        filename = directory + '/' + filename

    if not os.path.isfile(filename):
        py_printf('ERROR', 'Unable to parse convergence data since "{0}" is ' +
                  'not an existing OpenMOC log file'.format(filename))

    # Compile regular expressions to find the residual and eigenvalue data
    res = re.compile(b'res = ([0-9].[0-9]+E[+|-][0-9]+)')
    keff = re.compile(b'k_eff = ([0-9]+.[0-9]+)')

    # Parse the eigenvalues
    with open(filename, 'r+') as f:
        data = mmap.mmap(f.fileno(), 0)
        eigenvalues = keff.findall(data)

    # Parse the source residuals
    with open(filename, 'r+') as f:
        data = mmap.mmap(f.fileno(), 0)
        residuals = res.findall(data)

    # Create NumPy arrays of the data
    eigenvalues = np.array([float(eigenvalue) for eigenvalue in eigenvalues])
    residuals = np.array([float(residual) for residual in residuals])

    # Find the total number of source iterations
    num_iters = len(residuals)

    # Store the data in a dictionary to return to the user
    convergence_data = dict()
    convergence_data['# iters'] = num_iters
    convergence_data['eigenvalues'] = eigenvalues
    convergence_data['residuals'] = residuals
    return convergence_data


class Mesh(object):
    """A structured Cartesian mesh in two or three dimensions

    Attributes
    ----------
    dimension : Iterable of int
        The number of mesh cells in each direction.
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinate are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    width : Iterable of float
        The width of mesh cells in each direction.
    num_mesh_cells : Integral
        The total number of mesh cells in all dimensions
    mesh_cell_volume : Real
        The volume or area of a mesh cell in 2D or 3D in units of cm^3 or cm^2

    """

    def __init__(self):
        self._dimension = None
        self._lower_left = None
        self._upper_right = None
        self._width = None

    @property
    def dimension(self):
        return self._dimension

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @property
    def width(self):
        return self._width

    @property
    def num_mesh_cells(self):
        return np.prod(self._dimension)

    @property
    def mesh_cell_volume(self):
        return reduce(operator.mul, self.dimension, 1)

    @dimension.setter
    def dimension(self, dimension):
        cv.check_type('mesh dimension', dimension, Iterable, Integral)
        cv.check_length('mesh dimension', dimension, 2, 3)
        self._dimension = dimension

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('mesh lower_left', lower_left, Iterable, Real)
        cv.check_length('mesh lower_left', lower_left, 2, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        cv.check_type('mesh upper_right', upper_right, Iterable, Real)
        cv.check_length('mesh upper_right', upper_right, 2, 3)
        self._upper_right = upper_right

    @width.setter
    def width(self, width):
        cv.check_type('mesh width', width, Iterable, Real)
        cv.check_length('mesh width', width, 2, 3)
        self._width = width


    def get_mesh_cell_indices(self, point):
        """Get the mesh cell indices for a point within the geometry.

        Parameters
        ----------
        point : openmoc.Point
            A point within the geometry

        Returns
        -------
        indices : 2- or 3-tuple of Integral
            The mesh cell indices for the point. If the mesh is 2D then indices
            for x and y are returned; if the mesh is 3d indices for x, y, and z
            are returned.

        """

        cv.check_type('point', point, openmoc.Point)

        # Extract the x,y,z coordinates from the OpenMOC Point
        x, y, z = point.getX(), point.getY(), point.getZ()

        # Translate the point with respect to the center of the mesh
        x -= (self.upper_right[0] + self.lower_left[0]) / 2.
        y -= (self.upper_right[1] + self.lower_left[1]) / 2.
        if len(self.dimension) != 2:
            z -= (self.upper_right[2] + self.lower_left[2]) / 2.

        # Compute the mesh cell indices
        mesh_x = (x + self.dimension[0] * self.width[0] * 0.5) / self.width[0]
        mesh_y = (y + self.dimension[1] * self.width[1] * 0.5) / self.width[1]
        if len(self.dimension) == 2:
            mesh_z = 0
        else:
            mesh_z = (z + self.dimension[2] * self.width[2] * 0.5) / self.width[2]

        # Round the mesh cell indices down
        mesh_x = int(math.floor(mesh_x))
        mesh_y = int(math.floor(mesh_y))
        mesh_z = int(math.floor(mesh_z))

        # Throw error if indices are outside of the Mesh
        if len(self.dimension) == 2:
            if (mesh_x < 0 or mesh_x >= self.dimension[0]) or \
               (mesh_y < 0 or mesh_y >= self.dimension[1]):
                return np.nan, np.nan, np.nan
        else:
            if (mesh_x < 0 or mesh_x >= self.dimension[0]) or \
               (mesh_y < 0 or mesh_y >= self.dimension[1]) or \
               (mesh_z < 0 or mesh_z >= self.dimension[2]):
                return np.nan, np.nan, np.nan

        # Return mesh cell indices
        if len(self.dimension) == 2:
            return mesh_x, mesh_y
        else:
            return mesh_x, mesh_y, mesh_z

    @classmethod
    def from_lattice(cls, lattice, division=1):
        """Create a mesh from an existing lattice

        Parameters
        ----------
        lattice : openmoc.Lattice
            Uniform rectangular lattice used as a template for this mesh.
        division : int
            Number of mesh cells per lattice cell.
            If not specified, there will be 1 mesh cell per lattice cell.

        Returns
        -------
        openmoc.process.Mesh
            Mesh instance

        """
        cv.check_type("lattice", lattice, openmoc.Lattice)
        cv.check_type("division", division, Integral)
        if lattice.getNonUniform():
           raise ValueError("Lattice must be uniform.")

        shape = np.array((lattice.getNumX(), lattice.getNumY(),
                          lattice.getNumZ()))
        width = np.array((lattice.getWidthX(), lattice.getWidthY(),
                          lattice.getWidthZ()))
        lleft = np.array((lattice.getMinX(), lattice.getMinY(),
                          lattice.getMinZ()))
        uright = lleft + shape*width
        uright[np.isinf(width)] = np.inf

        mesh = cls()
        mesh.width = width
        mesh.lower_left = lleft
        mesh.upper_right = uright
        mesh.dimension = [s*division for s in shape]

        return mesh

    def tally_fission_rates(self, solver, volume='integrated', nu=False):
        """Compute the fission rates in each mesh cell.

        NOTE: This method assumes that the mesh perfectly aligns with the
        flat source region mesh used in the OpenMOC calculation.

        NOTE: The user must supply 'fission' as well as 'nu-fission' multi-group
        cross sections to each material in the geometry. Although 'nu-fission'
        is all that is required for an MOC calculation, 'fission' is what is
        used to compute the fission rates.

        Parameters
        ----------
        solver : {openmoc.CPUSolver, openmoc.GPUSolver, openmoc.VectorizedSolver}
            The solver used to compute the flux
        volume : {'averaged' ,'integrated'}
            Compute volume-averaged or volume-integrated fission rates
        nu : bool
            Find 'nu-fission' rates instead of 'fission' rates

        Returns
        -------
        tally : numpy.ndarray of Real
            A NumPy array of the fission rates tallied in each mesh cell

        """

        global solver_types
        cv.check_type('solver', solver, solver_types)
        cv.check_value('volume', volume, ('averaged', 'integrated'))

        geometry = solver.getGeometry()
        num_fsrs = int(geometry.getNumTotalFSRs())

        # Compute the volume- and energy-integrated fission rates for each FSR
        fission_rates = \
            solver.computeFSRFissionRates(int(geometry.getNumTotalFSRs()), nu)

        # Initialize a 2D or 3D NumPy array in which to tally
        tally = np.zeros(tuple(self.dimension), dtype=np.float)

        # Tally the fission rates in each FSR to the corresponding mesh cell
        for fsr in range(num_fsrs):
            point = geometry.getFSRPoint(fsr)
            mesh_indices = self.get_mesh_cell_indices(point)

            if np.nan not in mesh_indices:
                tally[mesh_indices] += fission_rates[fsr]

        # Average the fission rates by mesh cell volume if needed
        if volume == 'averaged':
            tally /= self.mesh_cell_volume

        return tally


    def tally_reaction_rates_on_mesh(self, solver, rxn_type,
                                     volume='integrated', energy='integrated'):
        """Compute 'material' or 'cell' reaction rates on a mesh

        This method streamlines the process of tallying reaction rates on a
        mesh by constructing the `domains_to_coeffs' dictionary and wrapping
        the Mesh.tally_on_mesh() method.

        NOTE: This method assumes that the mesh perfectly aligns with the
        flat source region mesh used in the OpenMOC calculation.

        Parameters
        ----------
        solver : {openmoc.CPUSolver, openmoc.GPUSolver, openmoc.VectorizedSolver}
            The solver used to compute the flux
        rxn_type : {'flux', 'total', 'scatter',
                'fission', 'nu-fission'}
        volume : {'averaged', 'integrated'}
            Compute volume-averaged or volume-integrated tallies
        energy : {'by_group', 'integrated'}
            Compute tallies by energy group or integrate across groups

        Returns
        -------
        tally : numpy.ndarray of Real
            A NumPy array of the fission rates tallied in each mesh cell indexed
            by FSR ID and energy group (if energy is 'by_group')

        """
        global rxn_types
        cv.check_value('rxn_type', rxn_type, rxn_types)

        geometry = solver.getGeometry()
        num_groups = geometry.getNumEnergyGroups()
        # Functionally, "material" and "cell" will produce identical results.
        # It only affects how we build the `domains_to_coeffs' dictionary.
        domain_type = "material"
        domains_to_coeffs = {}
        for k, mat in geometry.getAllMaterials().items():
            domains_to_coeffs[k] = np.zeros(num_groups)
            for g in range(num_groups):
                sigma = get_sigma_by_group(mat, rxn_type, g)
                domains_to_coeffs[k][g] = sigma
        tally = self.tally_on_mesh(solver, domains_to_coeffs, domain_type,
                                   volume, energy)
        return tally


    def tally_on_mesh(self, solver, domains_to_coeffs, domain_type='fsr',
                      volume='integrated', energy='integrated'):
        """Compute arbitrary reaction rates in each mesh cell.

        NOTE: This method assumes that the mesh perfectly aligns with the
        flat source region mesh used in the OpenMOC calculation.

        Parameters
        ----------
        solver : {openmoc.CPUSolver, openmoc.GPUSolver, openmoc.VectorizedSolver}
            The solver used to compute the flux
        domains_to_coeffs : dict or numpy.ndarray of Real
            A mapping of spatial domains and energy groups to the coefficients
            to multiply the flux in each domain. If domain_type is 'material'
            or 'cell' then the coefficients must be a Python dictionary indexed
            by material/cell ID mapped to NumPy arrays indexed by energy group.
            If domain_type is 'fsr' then the coefficients may be a dictionary
            or NumPy array indexed by FSR ID and energy group. Note that the
            energy group indexing should start at 0 rather than 1 for the
            highest energy in accordance with Python's 0-based indexing.
        domain_type : {'fsr', 'cell', 'material'}
            The type of domain for which the coefficients are defined
        volume : {'averaged', 'integrated'}
            Compute volume-averaged or volume-integrated tallies
        energy : {'by_group', 'integrated'}
            Compute tallies by energy group or integrate across groups

        Returns
        -------
        tally : numpy.ndarray of Real
            A NumPy array of the fission rates tallied in each mesh cell indexed
            by FSR ID and energy group (if energy is 'by_group')

        """

        global solver_types
        cv.check_type('solver', solver, solver_types)
        cv.check_value('domain_type', domain_type, ('fsr', 'cell', 'material'))
        cv.check_value('volume', volume, ('averaged', 'integrated'))
        cv.check_value('energy', energy, ('by_group', 'integrated'))

        # Extract parameters from the Geometry
        geometry = solver.getGeometry()
        num_groups = geometry.getNumEnergyGroups()
        num_fsrs = int(geometry.getNumTotalFSRs())

        # Coefficients must be specified as a dict, ndarray or DataFrame
        if domain_type in ['material', 'cell']:
            cv.check_type('domains_to_coeffs', domains_to_coeffs, dict)
        else:
            cv.check_type('domains_to_coeffs',
                          domains_to_coeffs, (dict, np.ndarray))

        # Extract the FSR fluxes from the Solver
        fluxes = get_scalar_fluxes(solver)

        # Initialize a 2D or 3D NumPy array in which to tally
        tally_shape = tuple(self.dimension) + (num_groups,)
        tally = np.zeros(tally_shape, dtype=np.float)

        # Compute product of fluxes with domains-to-coeffs mapping by group, FSR
        for fsr in range(num_fsrs):
            point = geometry.getFSRPoint(fsr)
            mesh_indices = self.get_mesh_cell_indices(point)

            if np.nan in mesh_indices:
                continue

            volume = solver.getFSRVolume(fsr)
            fsr_tally = np.zeros(num_groups, dtype=np.float)

            # Determine domain ID (material, cell or FSR) for this FSR
            if domain_type == 'fsr':
                domain_id = fsr
            else:
                coords = \
                    openmoc.LocalCoords(point.getX(), point.getY(), point.getZ())
                coords.setUniverse(geometry.getRootUniverse())
                cell = geometry.findCellContainingCoords(coords)
                if domain_type == 'cell':
                    domain_id = cell.getId()
                else:
                    domain_id = cell.getFillMaterial().getId()

            # Tally flux multiplied by coefficients by energy group
            for group in range(num_groups):
                fsr_tally[group] = \
                    fluxes[fsr, group] * domains_to_coeffs[domain_id][group]

            # Increment mesh tally with volume-integrated FSR tally
            tally[mesh_indices] += fsr_tally * volume

        # Integrate the energy groups if needed
        if energy == 'integrated':
            tally = np.sum(tally, axis=len(self.dimension))

        # Average the fission rates by mesh cell volume if needed
        if volume == 'averaged':
            tally /= self.mesh_cell_volume

        return tally
