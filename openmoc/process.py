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
#         openmoc.materialize module.
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
import os

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from log import *
# For Python 3.X.X
else:
  from openmoc.log import *


##
# @brief This routine is computes the fission rate in each flat source region,
#        and combines the rates into pin-wise and assembly-wise fission rates
#        in an exported binary file.
# @details This routine is intended to be called by the user in Python to
#          compute pin and assembly fission rates. The routine either exports
#          fission rates to an HDF5 binary file or ASCII files for each nested
#          Universe hierarchy. The routine makes use of a recursive utility
#          function computeUniverseFissionRate(...) which is also implemented
#          in the process module. This routine may be called from a Python
#          script as follows:
#
# @code
#          compute_pin_powers(solver, use_hdf5=True)
# @endcode
#
# @param solver a pointer to a Solver class object
# @param use_hdf5 whether or not to export pin powers to an HDF5 file
def compute_pin_powers(solver, use_hdf5=False):

  # Error checking of input parameters

  directory = get_output_directory() + '/pin-powers/'

  # Determine which Universes and Lattices contain fissionable Materials
  geometry = solver.getGeometry()
  geometry.computeFissionability()

  # Compute the volume-weighted fission rates for each FSR
  fission_rates = solver.computeFSRFissionRates(geometry.getNumFSRs())

  # Initialize a new HDF5 file to store the pin power data
  if use_hdf5:

    import h5py

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    f = h5py.File(directory + 'fission-rates.h5', 'w')
    f.close()

  # Get the base Universe in the Geometry and compute pin powers for each
  # level of nested Universes and Lattices
  universe = geometry.getUniverse(0)
  compute_universe_fission_rate(geometry, universe, 0, fission_rates,
                                use_hdf5, directory=directory)


##
# @brief A recursive routine to compute the fission rate for all cells in a
#        given universe, and for all universes contained within those cells.
# @details This method is an internal utility helper method for the
#          compute_pin_powers(...) routine and is NOT intended to be
#          called by the user.
# @param geometry a pointer to a Geometry object
# @param universe a pointer to the universe of interest
# @param FSR_offset the offset for this universe in the nested hierarchy
# @param FSR_fission_rates the array of fission rates for all FSRs
# @param use_hdf5 whether or not to export pin powers to HDF5
# @param attributes a list of strings for the path in the nested CSG hierarchy
# @param directory the output directory for this universe's fission rate
def compute_universe_fission_rate(geometry, universe, FSR_offset,
                                  FSR_fission_rates, use_hdf5=False,
                                  attributes = [], directory=''):

  fission_rate = 0.0

  # If the Universe is not fissionable, the fission rate within it is zero
  if not universe.isFissionable():
    return fission_rate

  # If the Universe is a fissionable SIMPLE type Universe
  elif universe.getType() is SIMPLE:

    # Create a directory/file for this Universe's total fission rate
    attributes.append('universe' + str(universe.getId()))

    num_cells = universe.getNumCells()
    cell_ids = universe.getCellIds(int(num_cells))

    # For each of the Cells inside the Universe, check if it is a
    # MATERIAL or FILL type
    for cell_id in cell_ids:

      cell = universe.getCell(int(cell_id))

      # If the current cell is a MATERIAL type cell
      if cell.getType() is MATERIAL:
        cell = castCellToCellBasic(cell)
        fsr_id = universe.getFSR(cell.getId()) + FSR_offset
        fission_rate += FSR_fission_rates[fsr_id]

      # The current Cell is a FILL type cell
      else:
        cell = castCellToCellFill(cell)
        universe_fill = cell.getUniverseFill()
        fsr_id = universe.getFSR(cell.getId()) + FSR_offset
        fission_rate += \
          compute_universe_fission_rate(geometry, universe_fill, fsr_id, \
                                        FSR_fission_rates, use_hdf5, \
                                        attributes, directory)

    # Remove the subdirectory structure attribute for this Universe
    attributes.pop()


  # This is a fissionable LATTICE type Universe
  else:

    # Add a new attribute for this Lattice in the pin power hierarchy
    attributes.append('lattice' + str(universe.getId()))

    # Retrieve a pointer to this Lattice and its dimensions
    lattice = castUniverseToLattice(universe)
    num_x = lattice.getNumX()
    num_y = lattice.getNumY()

    # Create a NumPy array to store all of this Lattice's pin powers
    lattice_cell_powers = np.zeros((num_x, num_y))

    attributes.append('x')
    attributes.append('y')

    # Loop over all Lattice cells in this Lattice
    for i in range(num_y-1, -1, -1):
      for j in range(num_x):

        # Get a pointer to the current Lattice cell
        cell_universe = lattice.getUniverse(j,i)

        # Get the FSR Id prefix for this Lattice cell
        fsr_id = lattice.getFSR(j,i) + FSR_offset

        attributes[-1] = 'y' + str(i)
        attributes[-2] = 'x' + str(j)

        # Compute the fission rate within this Lattice cell
        lattice_cell_powers[j,i] = \
            compute_universe_fission_rate(geometry, cell_universe, fsr_id,
                                          FSR_fission_rates, use_hdf5,
                                          attributes, directory)

        # Increment total Lattice power
        fission_rate += lattice_cell_powers[j,i]

    # Remove subdirectory attributes for x, y Lattice cell indices
    attributes.pop()
    attributes.pop()
    attributes.pop()

    # Flip the x-dimension of the array so that it is indexed
    # starting from the top left corner as
    lattice_cell_powers = np.fliplr(lattice_cell_powers)


    #######################################################################
    #######################  STORE PIN CELL DATA  #########################
    #######################################################################

    # If using HDF5, store data categorized by attributes in an HDF5 file
    if use_hdf5:

      import h5py

      # Create a h5py file handle for the file
      f = h5py.File(directory + 'fission-rates.h5', 'a')

      # Create the group for this Lattice, Universe combination
      curr_group = f.require_group(str.join('/', attributes))

      # Create a new dataset for this Lattice's pin powers
      curr_group.create_dataset('fission-rates', data=lattice_cell_powers)

      f.close()


    # If not using HDF5, store data categorized by subdirectory in ASCII files
    else:

      subdirectory = directory + str.join('/', attributes) + '/'

      # Make directory if it does not exist
      if not use_hdf5 and not os.path.exists(subdirectory):
        os.makedirs(subdirectory)

        np.savetxt(subdirectory + 'fission-rates.txt',
                   lattice_cell_powers, delimiter=',')

  # Return the total fission rate for this Lattice
  return fission_rate



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
#          store_simulation_state(solver, fluxes=True, sources=True, \
#                                 pin_powers=True, use_hdf5=True)
# @endcode
#
# @param solver a pointer to a Solver object
# @param fluxes whether to store FSR scalar fluxes (false by default)
# @param sources whether to store FSR sources (false by default)
# @param pin_powers whether to store fission rates (false by default)
# @param use_hdf5 whether to export to HDF5 (default) or Python pickle file
# @param filename the filename to use (default is 'simulation-state.h5')
# @param append append to existing file or create new one (false by default)
# @param note an additional string note to include in state file
def store_simulation_state(solver, fluxes=False, sources=False,
                           pin_powers=False, use_hdf5=False,
                           filename='simulation-state',
                           append=True, note=''):

  import datetime

  directory = 'simulation-states'

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
    sources = np.zeros((num_FSRs, num_groups))

    # Get the scalar flux for each FSR and energy group
    for i in range(num_FSRs):
      for j in range(num_groups):
        sources[i,j] = solver.getFSRSource(i,j+1)

  # If the user requested to store pin powers
  if pin_powers:

    # Generate and store pin powers
    compute_pin_powers(solver, use_hdf5=use_hdf5)

  # If using HDF5
  if use_hdf5:

    import h5py

    # Create a file handle
    if append:
      f = h5py.File(directory + '/' + filename + '.h5', 'a')
    else:
      f = h5py.File(directory + '/' + filename + '.h5', 'w')

    # Create groups for the day and time in the HDF5 file
    day_group = f.require_group(str(month)+'-'+str(day)+'-'+str(year))
    time_group = day_group.create_group(str(hr)+':'+str(mins)+':'+str(sec))

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
      time_group.create_dataset('FSR sources', data=sources)

    if pin_powers:

      # Open the pin powers file generated by compute_pin_powers(...)
      pin_powers_file = h5py.File('pin-powers/fission-rates.h5', 'r')

      # Deep copy the group of pin powers from pin_powers_file
      f.copy(pin_powers_file, time_group, name='fission-rates')

      # Close the pin powers file
      pin_powers_file.close()

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
    day = str(month)+'-'+str(day)+'-'+str(year)
    time = str(hr)+':'+str(mins)+':'+str(sec)

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
      state['FSR sources'] = sources

    if powers:
      py_printf('WARNING', 'The process.storeSimulationState(...)' + \
                'method only supports pin power storage for HDF5 files')

    # Pickle the simulation states to a file
    pickle.dump(sim_states, open(filename, 'wb'))


##
# @brief
# @details

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
#          Note: If the pin and and assembly fission rates were stored to
#          the binary file, they are not restored and returned in this method.
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

  # If using HDF5
  if '.h5' in filename or '.hdf5' in filename:

    import h5py

    # Create a file handle
    f = h5py.File(directory + '/' + filename, 'r')

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
          sources = dataset['FSR sources'][...]
          state['FSR sources'] = sources

        if 'note' in dataset:
          state['note'] = str(dataset['note'])

        if 'fission-rates' in dataset:
          py_printf('WARNING', 'The process.restore_simulation_state(...)'+\
                    'method does not yet support pin powers')

    return states


  # If using a Python pickled file
  elif '.pkl' in filename:

    import pickle

    # Load the dictionary from the pickle file
    filename = directory + '/' + filename
    states = pickle.load(file(filename, 'rb'))

    return states


  # If file does not have a recognizable extension
  else:
    py_printf('WARNING', 'Unable to restore the simulation states file %s' + \
              ' since it does not have a supported file extension. Only ' + \
              '*.h5, *.hdf5, and *.pkl files are supported', filename)

    return {}
