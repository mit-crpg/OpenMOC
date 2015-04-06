##
# @file process.py
# @package openmoc.process
# @brief The process module provides utility functions for processing of
#        OpenMOC simulations.
# @details This module provides downstream data processing capabilities for
#          OpenMOC simulations, such as storing and retrieving simulation
#          statepoint files, computing pin/assembly powers, and more.
# @author William Boyd (wboyd@mit.edu)
# @date April 27, 2013

import sys

## @var openmoc
#  @brief The openmoc module in use in the Python script using the
#         openmoc.process module.
openmoc = ''

# Determine which OpenMOC module is being used
if 'openmoc.gnu.double' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.double']
elif 'openmoc.gnu.single' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.single']
elif 'openmoc.intel.double' in sys.modules:
  openmoc = sys.modules['openmoc.intel.double']
elif 'openmoc.intel.single' in sys.modules:
  openmoc = sys.modules['openmoc.intel.single']
elif 'openmoc.bgq.double' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.double']
elif 'openmoc.bgq.single' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.single']
else:
  from openmoc import *

import numpy as np
import os, re, mmap

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from log import *
# For Python 3.X.X
else:
  from openmoc.log import *


##
# @brief This routine checks if a given value is an integer data type.
#
# @param val a value to check
def is_integer(val):
  return isinstance(val, (int, np.int32, np.int64))


##
# @brief This routine checks if a given value is a string data type.
#
# @param val a value to check
def is_string(val):
  return isinstance(val, (str, np.str))


##
# @brief This routine checks if a given value is an float data type.
#
# @param val a value to check
def is_float(val):
  return isinstance(val, (float, np.float32, np.float64))


##
# @brief This routine computes the fission rate in each flat source region,
#        and combines the rates based on their hierarchical universe/lattice
#        structure. The fission rates are then exported to a binary HDF5
#        or python pickle file.
# @details This routine is intended to be called by the user in Python to
#          compute fission rates. Typically, the fission rates will represent
#          pin powers. The routine either exports fission rates to an HDF5
#          binary file or pickle file with each fission rate being indexed by
#          a string representing the universe/lattice hierarchy.
#          This routine may be called from a Python script as follows:
#
# @code
#          compute_fission_rates(solver, use_hdf5=True)
# @endcode
#
# @param solver a pointer to a Solver class object
# @param use_hdf5 whether or not to export fission rates to an HDF5 file
def compute_fission_rates(solver, use_hdf5=False):

  # create directory and filename
  directory = get_output_directory() + '/fission-rates/'
  filename = 'fission-rates'

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Get geometry
  geometry = solver.getGeometry()

  # Compute the volume-weighted fission rates for each FSR
  fsr_fission_rates = solver.computeFSRFissionRates(geometry.getNumFSRs())

  # Initialize fission rates dictionary
  fission_rates_sum = {}

  # Loop over FSRs and populate fission rates dictionary
  for fsr in range(geometry.getNumFSRs()):

    if geometry.findFSRMaterial(fsr).isFissionable():

      # Get the linked list of LocalCoords
      point = geometry.getFSRPoint(fsr)
      coords = LocalCoords(point.getX(), point.getY())
      coords.setUniverse(geometry.getRootUniverse())
      geometry.findCellContainingCoords(coords)
      coords = coords.getHighestLevel().getNext()

      # initialize dictionary key
      key = 'UNIV = 0 : '

      # Parse through the linked list and create fsr key.
      # If lowest level sub dictionary already exists, then increment 
      # fission rate; otherwise, set the fission rate.
      while True:
        if coords.getType() is LAT:
          key += 'LAT = ' + str(coords.getLattice().getId()) + ' (' + \
                 str(coords.getLatticeX()) + ', ' + \
                 str(coords.getLatticeY()) + ') : '
        else:
          key += 'UNIV = ' + str(coords.getUniverse().getId()) + ' : '

        # Remove the trailing ' : ' on the end of the key if at last univ/lat
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

  # If using HDF5
  if use_hdf5:

    import h5py

    # Open HDF5 file
    f = h5py.File(directory + filename + '.h5', 'w')

    # Write the fission rates to the HDF5 file
    fission_rates_group = f.create_group('fission-rates')
    for key, value in fission_rates_sum.items():
      fission_rates_group.attrs[key] = value

    # Close hdf5 file
    f.close()

  else:

    import pickle

    # Pickle the fission rates to a file
    pickle.dump(fission_rates_sum, open(directory + filename + '.pkl', 'wb'))

##
# @brief This method stores all of the data for an OpenMOC simulation to a
#        a binary file for downstream data processing.
# @details The method may be used to store the type of Solver used, floating
#          point precision, exponential evaluation method, number of FSRs,
#          number of materials, number of energy groups, number of azimuthal
#          angles, number of polar angles, track spacing, number of Tracks,
#          number of Track segments, number of source iterations, source
#          convergence tolerance, converged \f$ k_{eff} \f$, total runtime,
#          and number of OpenMP or CUDA threads. In addition, the routine
#          can store the FSR flux array, FSR source array, and pin and
#          assembly fission rates.
#
#          The routine may export the simulation data to either an HDF5 or
#          a Python pickle binary file. Users may tell the routine to either
#          create a new binary output file, or append to an existing file
#          using a timestamp to record multiple simulation states to the
#          same file.
#
#          This method may be called from Python as follows:
#
# @code
#          store_simulation_state(solver, fluxes=True, source=True, \
#                                 fission_rates=True, use_hdf5=True)
# @endcode
#
# @param solver a pointer to a Solver object
# @param fluxes whether to store FSR scalar fluxes (false by default)
# @param sources whether to store FSR sources (false by default)
# @param fission_rates whether to store fission rates (false by default)
# @param use_hdf5 whether to export to HDF5 (default) or Python pickle file
# @param filename the filename to use (default is 'simulation-state.h5')
# @param directory the directory to use (default is 'simulation-states')
# @param append append to existing file or create new one (false by default)
# @param note an additional string note to include in state file
def store_simulation_state(solver, fluxes=False, sources=False,
                           fission_rates=False, use_hdf5=False,
                           filename='simulation-state', 
                           directory = 'simulation-states',
                           append=True, note=''):

  import datetime

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
  elif 'ThreadPrivateSolver' in str(solver.__class__):
    solver_type = 'ThreadPrivateSolver'
  elif 'VectorizedSolver' in str(solver.__class__):
    solver_type = 'VectorizedSolver'
  elif 'VectorizedPrivateSolver' in str(solver.__class__):
    solver_type = 'VectorizedPrivateSolver'
  elif 'GPUSolver' in str(solver.__class__):
    solver_type = 'GPUSolver'

  # Determine the floating point precision level
  if solver.isUsingDoublePrecision():
    precision = 'double'
  else:
    precision = 'single'

  # Determine whether we are using the exponential intrinsic or
  # linear interpolation for exponential evaluations
  if solver.isUsingExponentialIntrinsic():
      method = 'exp intrinsic'
  else:
    method = 'linear interpolation'

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
  num_FSRs = geometry.getNumFSRs()
  num_materials = geometry.getNumMaterials()
  num_groups = geometry.getNumEnergyGroups()
  num_tracks = track_generator.getNumTracks()
  num_segments = track_generator.getNumSegments()
  spacing = track_generator.getTrackSpacing()
  num_azim = track_generator.getNumAzim()
  num_polar = solver.getNumPolarAngles()
  num_iters = solver.getNumIterations()
  thresh = solver.getSourceConvergenceThreshold()
  tot_time = solver.getTotalTime()
  keff = solver.getKeff()

  if solver_type is 'GPUSolver':
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
        scalar_fluxes[i,j] = solver.getFSRScalarFlux(i,j+1)

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

    import h5py

    # Create a file handle
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
    if not note is '':
      time_group.attrs['note'] = note

    # Store simulation data to the HDF5 file
    time_group.create_dataset('solver type', data=solver_type)
    time_group.create_dataset('# FSRs', data=num_FSRs)
    time_group.create_dataset('# materials', data=num_materials)
    time_group.create_dataset('# energy groups', data=num_groups)
    time_group.create_dataset('# tracks', data=num_tracks)
    time_group.create_dataset('# segments', data=num_segments)
    time_group.create_dataset('track spacing [cm]', data=spacing)
    time_group.create_dataset('# azimuthal angles', data=num_azim)
    time_group.create_dataset('# polar angles', data=num_polar)
    time_group.create_dataset('# iterations', data=num_iters)
    time_group.create_dataset('source residual threshold', data=thresh)
    time_group.create_dataset('exponential', data=method)
    time_group.create_dataset('floating point', data=precision)
    time_group.create_dataset('CMFD', data=cmfd)
    time_group.create_dataset('time [sec]', data=tot_time)
    time_group.create_dataset('keff', data=keff)

    if solver_type is 'GPUSolver':
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

      # Open the fission rates file generated by compute_fission_rates(...)
      fission_rates_file = h5py.File('fission-rates/fission-rates.h5', 'r')

      # Deep copy the group of fission rates from fission_rates_file
      f.copy(fission_rates_file, time_group, name='fission-rates')

      # Close the pin powers file
      fission_rates_file.close()

    # Close the HDF5 file
    f.close()

  # If not using HDF5, we are pickling all of the data
  else:

    import pickle

    # Load the dictionary from the Pickle file
    filename = directory + '/' + filename + '.pkl'
    if os.path.exists(filename) and append:
      sim_states = pickle.load(file(filename, 'rb'))
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
    if not note is '':
      state['note'] = note

    # Store simulation data to a Python dictionary
    state['solver type'] = solver_type
    state['# FSRs'] = num_FSRs
    state['# materials'] = num_materials
    state['# energy groups'] = num_groups
    state['# tracks'] = num_tracks
    state['# segments'] = num_segments
    state['track spacing [cm]'] = spacing
    state['# azimuthal angles'] = num_azim
    state['# polar angles'] = num_polar
    state['# iterations'] = num_iters
    state['source residual threshold'] = thresh
    state['exponential'] = method
    state['floating point'] = precision
    state['CMFD'] = cmfd
    state['time [sec]'] = tot_time
    state['keff'] = keff

    if solver_type is 'GPUSolver':
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
        pickle.load(file('fission-rates/fission-rates.pkl', 'rb'))

    # Pickle the simulation states to a file
    pickle.dump(sim_states, open(filename, 'wb'))

    # Pickle the simulation states to a file
    pickle.dump(sim_states, open(filename, 'wb'))


##
# @brief This method restores all of the data for an OpenMOC simulation from a
#        a binary file for downstream data processing to a Python dictionary.
# @details The routine may import the simulation state from either an HDF5 or a
#          Python pickle binary file created by the store_simulation_state(...)
#          method. The method may be used to restore the type of Solver used,
#          floating point precision, exponential evaluation method, number of
#          FSRs, number of materials, number of energy groups, number of
#          azimuthal angles, number of polar angles, track spacing, number of
#          Tracks, number of Track segments, number of source iterations, source
#          convergence tolerance, converged \f$ k_{eff} \f$, total runtime,
#          and number of OpenMP or CUDA threads. In addition, the routine
#          can restore the FSR flux array, FSR source array.
#
#          Note: If the fission rates were stored in a hdf5 binary file, 
#          they are not restored and returned in this method.
#
#          This method may be called from Python as follows:
#
# @code
#          restore_simulation_state(filename='simulation-state-v1.3.h5')
# @endcode
#
# @param filename the simulation state filename string
# @param directory the directory where to find the simulation state file
# @return a Python dictionary of key/value pairs for simulation state data
def restore_simulation_state(filename='simulation-state.h5',
                             directory='simulation-states'):

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

        # Extract simulation data and store it in the sub-dictionary
        solver_type = str(dataset['solver type'])
        state['solver type'] = solver_type

        num_FSRs = int(dataset['# FSRs'][...])
        state['# FSRs'] = num_FSRs

        num_materials = int(dataset['# materials'][...])
        state['# materials'] = num_materials

        num_tracks = int(dataset['# tracks'][...])
        state['# tracks'] = num_materials

        num_segments = int(dataset['# segments'][...])
        state['# segments'] = num_segments

        spacing = int(dataset['track spacing [cm]'][...])
        state['track spacing [cm]'] = spacing

        num_azim = int(dataset['# azimuthal angles'][...])
        state['# azimuthal angles'] = num_azim

        num_polar = int(dataset['# polar angles'][...])
        state['# polar angles'] = num_polar

        num_iters = int(dataset['# iterations'][...])
        state['# iterations'] = num_iters

        thresh = float(dataset['source residual threshold'][...])
        state['source residual threshold'] = thresh

        method = str(dataset['exponential'][...])
        state['exponential'] = method

        precision = str(dataset['floating point'][...])
        state['floating point'] = precision

        cmfd = str(dataset['CMFD'][...])
        state['CMFD'] = cmfd

        time = float(dataset['time [sec]'][...])
        state['time [sec]'] = time

        keff =  float(dataset['keff'][...])
        state['keff'] = keff

        if solver_type is 'GPUSolver':
          num_threads = int(dataset['# threads per block'])
          num_blocks = int(dataset['# thread blocks'])
        else:
          num_threads = int(dataset['# threads'][...])

        if 'FSR scalar fluxes' in dataset:
          fluxes = dataset['FSR scalar fluxes'][...]
          state['FSR scalar fluxes'] = fluxes

        if 'FSR sources' in dataset:
          sources_array = dataset['FSR sources'][...]
          state['FSR sources'] = sources_array

        if 'note' in dataset:
          state['note'] = str(dataset['note'])

        if 'fission-rates' in dataset:
          py_printf('WARNING', 'The process.restore_simulation_state(...)'+\
                    'method does not yet support fission rates')

    return states


  # If using a Python pickled file
  elif '.pkl' in filename:

    import pickle

    # Load the dictionary from the pickle file
    states = pickle.load(file(filename, 'rb'))

    return states

  # If file does not have a recognizable extension
  else:
    py_printf('WARNING', 'Unable to restore the simulation states file %s' + \
              ' since it does not have a supported file extension. Only ' + \
              '*.h5, *.hdf5, and *.pkl files are supported', filename)

    return {}


##
# @brief Parse an OpenMOC log file to obtain a simulation's convergence data.
# @details This method compiles the eigenvalue and source residuals from each
#          iteration of an OpenMOC simulation. This data is inserted into a 
#          Python dictionary under the key names 'eigenvalues' and 'residuals',
#          along with an integer `# iters`, and returned to the user.
#
#          This method may be called from Python as follows:
# @code
#          parse_convergence_data(filename='openmoc-XX-XX-XXXX--XX:XX:XX.log')
# @endcode
#
# @param filename the OpenMOC log filename string
# @param directory the directory where to find the log file
# @return a Python dictionary of key/value pairs for convergence data
def parse_convergence_data(filename, directory=''):

  # If the user specified a directory
  if len(directory) > 0:
    filename = directory + '/' + filename

  if not os.path.isfile(filename):
    py_printf('ERROR', 'Unable to parse convergence data since "{0}" is ' + \
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
