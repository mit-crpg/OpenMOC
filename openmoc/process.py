##
# @file process.py
# @package openmoc.process
# @brief
# @author William Boyd (wboyd@mit.edu)
# @date April 27, 2013

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from log import *

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



##
# @brief
# @param
# @param
#
def computeFSRPinPowers(solver, use_hdf5=False):

    # Error checking of input parameters

    directory = openmoc.getOutputDirectory() + '/pin-powers/'

    # Determine which universes and lattices contain fissionable materials
    geometry = solver.getGeometry()
    geometry.computeFissionability()

    # Compute the volume-weighted FSR fission rates for each FSR
    fission_rates = solver.computeFSRFissionRates(geometry.getNumFSRs())

    # Initialize a new HDF5 file to store the pin power data
    if use_hdf5:

        import h5py

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        f = h5py.File(directory + 'fission-rates.h5', 'w')
        f.close()


    # Get the base universe in the geometry and compute pin powers for each
    # level of nested universes and lattices
    universe = geometry.getUniverse(0)
    computeUniverseFissionRate(geometry, universe, 0, fission_rates, use_hdf5,
                               directory=directory)


##
# @brief
# @param
# @param
# @param
# @param
# @param
# @return
#
def computeUniverseFissionRate(geometry, universe, FSR_id, 
                               FSR_fission_rates, use_hdf5=False, 
                               attributes = [], directory=''):

    fission_rate = 0.0

    # If the universe is not fissionable, the fission rate within it
    # is zero
    if not universe.isFissionable():
        return fission_rate

    # If the universe is a fissionable SIMPLE type universe
    elif universe.getType() is openmoc.SIMPLE:

        # Create a directory/file for this universe's total fission rate
        attributes.append('universe' + str(universe.getId()))

        num_cells = universe.getNumCells()
        cell_ids = universe.getCellIds(int(num_cells))

        # For each of the cells inside the universe, check if it is 
        # MATERIAL or FILL type
        for cell_id in cell_ids:
        
            cell = universe.getCell(int(cell_id))
            
            # If the current cell is a MATERIAL type cell, 
            if cell.getType() is openmoc.MATERIAL:
                cell = openmoc.castCellToCellBasic(cell)
                fsr_id = universe.getFSR(cell.getId()) + FSR_id
                fission_rate += FSR_fission_rates[fsr_id]
            
            # The current cell is a FILL type cell
            else:
                cell = openmoc.castCellToCellFill(cell)
                universe_fill = cell.getUniverseFill()
                fsr_id = universe.getFSR(cell.getId()) + FSR_id
                fission_rate += \
                    computeUniverseFissionRate(geometry, universe_fill, 
                                               fsr_id, FSR_fission_rates, 
                                               use_hdf5, attributes, directory)

        # Remove the subdirectory structure attribute for this universe
        attributes.pop()


    # This is a fissionable LATTICE type universe
    else:

        # Add a new attribute for this lattice in the pin power hierarchy
        attributes.append('lattice' + str(universe.getId()))

        # Retrieve a pointer to this lattice and its dimensions
        lattice = openmoc.castUniverseToLattice(universe)
        num_x = lattice.getNumX()
        num_y = lattice.getNumY()

        # Create a numpy array to store all of this lattice's pin powers
        lattice_cell_powers = np.zeros((num_x, num_y))

        attributes.append('x')
        attributes.append('y')

        # Loop over all lattice cells in this lattice
        for i in range(num_y-1, -1, -1):
            for j in range(num_x):

                # Get a pointer to the current lattice cell
                cell_universe = lattice.getUniverse(j,i)
                
                # Get the FSR Id prefix for this lattice cell
                fsr_id = lattice.getFSR(j,i) + FSR_id

                attributes[-1] = 'y' + str(i)
                attributes[-2] = 'x' + str(j)

                # Compute the fission rate within this lattice cell
                lattice_cell_powers[j,i] = \
                    computeUniverseFissionRate(geometry, cell_universe, fsr_id, 
                                               FSR_fission_rates, use_hdf5,
                                               attributes, directory)

                # Increment total lattice power
                fission_rate += lattice_cell_powers[j,i]

        # Remove subdirectory attributes for x, y lattice cell indices
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
            f = h5py.File(directory + 'fission-rates.hdf5', 'a')

            # Create the group for this lattice, universe combination
            curr_group = f.require_group(str.join('/', attributes))

            # Create a new dataset for this lattice's pin powers
            curr_group.create_dataset('fission-rates', data=lattice_cell_powers)

            f.close()


        # If not using HDF5, store data categorized by subdirectory in txt files
        else:

            subdirectory = directory + str.join('/', attributes) + '/'

            # Make directory if it does not exist
            if not use_hdf5 and not os.path.exists(subdirectory):
                os.makedirs(subdirectory)

            np.savetxt(subdirectory + 'fission-rates.txt', 
                          lattice_cell_powers, delimiter=',')

    # Return the total fission rate for this lattice
    return fission_rate



##
# @brief
# @param solver
# @param fluxes
# @param filename
def storeSimulationState(solver, fluxes=False, sources=False, use_hdf5=False, 
                         filename='simulation-state', append=True, note=''):

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

    # Determine the solver type
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

    # Get the geometry and track_generator from the solver
    geometry = solver.getGeometry()
    track_generator = solver.getTrackGenerator()

    # Retrieve useful data from the solver, geometry and track generator
    num_FSRs = geometry.getNumFSRs()
    num_materials = geometry.getNumMaterials()
    num_groups = geometry.getNumEnergyGroups()
    num_tracks = track_generator.getNumTracks()
    num_segments = track_generator.getNumSegments()
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
        time_group.create_dataset('# azimuthal angles', data=num_azim)
        time_group.create_dataset('# polar angles', data=num_polar)
        time_group.create_dataset('# iterations', data=num_iters)
        time_group.create_dataset('source residual threshold', data=thresh)
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
        state['# azimuthal angles'] = num_azim
        state['# polar angles'] = num_polar
        state['# iterations'] = num_iters
        state['source residual threshold'] = thresh
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

        # Pickle the simulation states to a file
        pickle.dump(sim_states, open(filename, 'wb'))


        
def restoreSimulationStates(filename='simulation-state.pkl', 
                            directory='simulation-states'):

    # If using HDF5
    if '.h5' in filename or '.hdf5' in filename:

        import h5py

        # Create a file handle
        f = h5py.File(directory + '/' + filename, 'r')

        states = {}

        for day in f.keys():

            states[day] = {}

            for time in f[day]:

                dataset = f[day][time]
                states[day][time] = {}
                state = states[day][time]

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

                num_azim = int(dataset['# azimuthal angles'][...])
                state['# azimuthal angles'] = num_azim

                num_polar = int(dataset['# polar angles'][...])
                state['# polar angles'] = num_polar

                num_iters = int(dataset['# iterations'][...])
                state['# iterations'] = num_iters

                thresh = float(dataset['source residual threshold'][...])
                state['source residual threshold'] = thresh

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

        return states


    # If using a Pickled file
    elif '.pkl':

        import pickle

        # Load the dictionary from the Pickle file
        filename = directory + '/' + filename
        states = pickle.load(file(filename, 'rb'))
        
        return states

    # If file does not have a recognizable extension
    else:

        return {}
