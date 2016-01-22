from __future__ import print_function

import filecmp
import hashlib
from optparse import OptionParser
import os
import shutil
import sys

sys.path.insert(0, 'openmoc')
import openmoc
import openmoc.process


class TestHarness(object):
    """General class for running OpenMOC regression tests."""
    def __init__(self):

        # Define default simulation parameters
        self.num_threads = 1
        self.spacing = 0.1
        self.num_azim = 4
        self.max_iters = 10
        self.tolerance = 1E-5

        self.input_set = None
        self.track_generator = None
        self.solver = None
        self.calc_mode = 'eigenvalue'
        self.res_type = openmoc.FISSION_SOURCE

        self.parser = OptionParser()
        self.parser.add_option('--update', dest='update',
                               action='store_true', default=False)

        self._opts = None
        self._args = None

    def create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        self.track_generator = \
            openmoc.TrackGenerator(geometry, self.num_azim, self.spacing)

    def create_solver(self):
        """Instantiate a CPUSolver."""
        self.solver = openmoc.CPUSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)

    def generate_tracks(self):
        """Generate Tracks and segments."""
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        (self._opts, self._args) = self.parser.parse_args()
        if self._opts.update:
            self._update_results()
        else:
            self._execute_test()

    def _execute_test(self):
        """Build geometry, ray trace, run calculation, and verify results."""
        try:
            self.run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def _update_results(self):
        """Update the results_true using the current version of OpenMOC."""
        try:
            self.run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def run_openmoc(self):
        """Run an OpenMOC eigenvalue or fixed source calculation."""

        openmoc.log.set_log_level('NORMAL')

        self.input_set.create_materials()
        self.input_set.create_geometry()
        self.create_trackgenerator()
        self.generate_tracks()
        self.create_solver()

        if self.calc_mode == 'eigenvalue':
            self.solver.computeEigenvalue(self.max_iters, res_type=self.res_type)
        elif self.calc_mode == 'fixed source':
            self.solver.computeFlux(self.max_iters, res_type=self.res_type)
        elif self.calc_mode == 'subcritical':
            self.solver.computeSource(self.max_iters, res_type=self.res_type)
        else:
            msg = 'Unable to run OpenMOC in mode {}'.format(self.calc_mode)
            raise ValueError(msg)

    def _get_results(self, keff=True, fluxes=False, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Write out the number of iterations
        num_iters = self.solver.getNumIterations()
        outstr = '# Iterations: {0}\n'.format(num_iters)

        # Write out the eigenvalue
        if keff and self.calc_mode == 'eigenvalue':
            keff = self.solver.getKeff()
            outstr += 'keff: {0:12.6E}'.format(keff)

        if fluxes:
            # Get the fluxes for each FSR and energy group
            fluxes = openmoc.process.get_scalar_fluxes(self.solver)

            # Create a list of the floating point flux values
            fluxes = ['{0:12.6E}'.format(flux) for flux in fluxes.ravel()]

            # Add the fluxes to the output string
            outstr += 'fluxes:\n'
            outstr += '\n'.join(fluxes) + '\n'

        # Hash the results if necessary.
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _write_results(self, results_string):
        """Write the results to an ASCII file."""
        with open('results_test.dat', 'w') as fh:
            fh.write(results_string)

    def _overwrite_results(self):
        """Overwrite the results_true with the results_test."""
        shutil.copyfile('results_test.dat', 'results_true.dat')

    def _compare_results(self):
        """Make sure the current results agree with the _true standard."""
        compare = filecmp.cmp('results_test.dat', 'results_true.dat')
        if not compare:
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree.'

    def _cleanup(self):
        """Delete track, plot, etc. directories and test files."""

        # Create list of directories and/or files to remove
        outputs = [os.path.join(os.getcwd(), 'tracks')]
        outputs.append(os.path.join(os.getcwd(), 'plots'))
        outputs.append(os.path.join(os.getcwd(), 'simulation-state'))
        outputs.append(os.path.join(os.getcwd(), 'results_test.dat'))

        # Remove each file / directory if it exists
        for output in outputs:
            if os.path.isfile(output):
                os.remove(output)
            elif os.path.isdir(output):
                shutil.rmtree(output)


class HashedTestHarness(TestHarness):
    """Specialized TestHarness that hashes the results."""

    def _get_results(self):
        """Digest info in the results and return as a string."""
        return super(HashedTestHarness, self)._get_results(hash_output=True)