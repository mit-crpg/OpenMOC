#!/usr/bin/env python

import os
import sys
import hashlib
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc.process as process
import numpy as np


class SimulationStateTestHarness(TestHarness):
    """An eigenvalue calculation with storage and retrieval of the simulation
    state using the openmoc.process module."""

    def __init__(self):
        super(SimulationStateTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _run_openmoc(self):
        """Run an eigenvalue calculation and store the simulation state."""
        super(SimulationStateTestHarness, self)._run_openmoc()

        # Store the simulatiation state in HDF5 and pickled binary files
        process.store_simulation_state(self.solver, append=False, fluxes=True,
                                       sources=True, use_hdf5=True)
        process.store_simulation_state(self.solver, append=False, fluxes=True,
                                       sources=True, use_hdf5=False)

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=False):
        """Digest info in the simulation states and return hash as a string."""

        # Restore simulation state from pickled binary file
        pkl_state = \
            process.restore_simulation_state(filename='simulation-state.pkl')
        hdf5_state = \
            process.restore_simulation_state(filename='simulation-state.h5')

        # Add key, value pairs from pickled simulation state to ouptut string
        outstr = 'HDF5 Simulation State\n'
        for key1 in sorted(hdf5_state):
            for key2 in sorted(hdf5_state[key1]):
                for key3 in sorted(hdf5_state[key1][key2]):
                    value = hdf5_state[key1][key2][key3]

                    # Ignore runtime since it changes for each run
                    if 'time' in key3:
                        continue
                    # Ignore thread counts for sequential/parallel runs
                    if 'threads' in key3:
                        continue

                    # Create a list of the floating point flux values
                    if isinstance(value, np.ndarray):
                        outstr += '{0}: '.format(key3)
                        values = \
                            ['{0:12.6E}'.format(val) for val in value.ravel()]
                        outstr += '\n'.join(values) + '\n'
                    elif isinstance(value, float):
                        outstr += '{0}:\t{1:.10f}\n'.format(key3, value)
                    else:
                        outstr += '{0}:\t{1}\n'.format(key3, value)

        # Add key, value pairs from HDF5 simulation state to ouptut string
        outstr += 'Pickle Simulation State\n'
        for key1 in sorted(pkl_state):
            for key2 in sorted(pkl_state[key1]):
                for key3 in sorted(pkl_state[key1][key2]):
                    value = pkl_state[key1][key2][key3]

                    # Ignore runtime since it changes for each run
                    if 'time' in key3:
                        continue
                    # Ignore thread counts for sequential/parallel runs
                    if 'threads' in key3:
                        continue

                    # Create a list of the floating point flux values
                    if isinstance(value, np.ndarray):
                        outstr += '{0}: '.format(key3)
                        values = \
                            ['{0:12.6E}'.format(val) for val in value.ravel()]
                        outstr += '\n'.join(values) + '\n'
                    elif isinstance(value, float):
                        outstr += '{0}:\t{1:.10f}\n'.format(key3, value)
                    else:
                        outstr += '{0}:\t{1}\n'.format(key3, value)

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


if __name__ == '__main__':
    harness = SimulationStateTestHarness()
    harness.main()
