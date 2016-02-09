from __future__ import print_function

import filecmp
import hashlib
from optparse import OptionParser
import os
import shutil
import sys
from collections import OrderedDict

sys.path.insert(0, 'openmoc')
import openmoc
import openmoc.process


class TestHarness(object):
    """General class for running OpenMOC regression tests."""
    def __init__(self):

        openmoc.log.set_log_level('NORMAL')

        # Define default simulation parameters
        self.spacing = 0.1
        self.num_azim = 4
        self.max_iters = 500
        self.tolerance = 1E-5

        # Define threads based on OMP_NUM_THREADS env var set by run_tests.py
        if 'OMP_NUM_THREADS' in os.environ:
            self.num_threads = int(os.environ['OMP_NUM_THREADS'])
        else:
            self.num_threads = 1

        self.input_set = None
        self.track_generator = None
        self.solver = None
        self.solution_type = 'eigenvalue'
        self.calculation_mode = openmoc.FORWARD
        self.res_type = openmoc.FISSION_SOURCE

        self.parser = OptionParser()
        self.parser.add_option('--update', dest='update',
                               action='store_true', default=False)

        self._opts = None
        self._args = None


    def _create_geometry(self):
        """Initialize the materials and geometry in the InputSet."""
        self.input_set.create_materials()
        self.input_set.create_geometry()

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        self.track_generator = \
            openmoc.TrackGenerator(geometry, self.num_azim, self.spacing)

    def _create_solver(self):
        """Instantiate a CPUSolver."""
        self.solver = openmoc.CPUSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)

    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Always use 1 thread for FSR reproducibility
        self.track_generator.setNumThreads(1)
        self.track_generator.generateTracks()

    def _setup(self):
        """Build materials, geometry and perform ray tracing."""
        self._create_geometry()
        self._create_trackgenerator()
        self._generate_tracks()
        self._create_solver()


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
            self._setup()
            self._run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def _update_results(self):
        """Update the results_true using the current version of OpenMOC."""
        try:
            self._setup()
            self._run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue or fixed source calculation."""

        if self.solution_type == 'eigenvalue':
            self.solver.computeEigenvalue(self.max_iters, res_type=self.res_type,
                                          mode=self.calculation_mode)
        elif self.solution_type == 'flux':
            self.solver.computeFlux(self.max_iters)
        elif self.solution_type == 'source':
            self.solver.computeSource(self.max_iters, res_type=self.res_type,
                                      mode=self.calculation_mode)
        else:
            msg = 'Unable to run OpenMOC in mode {0}'.format(self.solution_type)
            raise ValueError(msg)

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the solver and return as a string."""

        outstr = ''

        # Write out the number of iterations
        if num_iters:
            num_iters = self.solver.getNumIterations()
            outstr += '# Iterations: {0}\n'.format(num_iters)

        # Write out the eigenvalue
        if keff and self.solution_type == 'eigenvalue':
            keff = self.solver.getKeff()
            outstr += 'keff: {0:12.5E}\n'.format(keff)

        if fluxes:
            # Get the fluxes for each FSR and energy group
            fluxes = openmoc.process.get_scalar_fluxes(self.solver)

            # Create a list of the floating point flux values
            fluxes = ['{0:12.6E}'.format(flux) for flux in fluxes.ravel()]

            # Add the fluxes to the output string
            outstr += 'fluxes:\n'
            outstr += '\n'.join(fluxes) + '\n'

        # Write out the number of FSRs
        if num_fsrs:
            num_fsrs = self.input_set.geometry.getNumFSRs()
            outstr += '# FSRs: {0}\n'.format(num_fsrs)

        # Write out the number of tracks
        if num_tracks:
            num_tracks = self.track_generator.getNumTracks()
            outstr += '# tracks: {0}\n'.format(num_tracks)

        # Write out the number of segments
        if num_segments:
            num_segments = self.track_generator.getNumSegments()
            outstr += '# segments: {0}\n'.format(num_segments)

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
        outputs.append(os.path.join(os.getcwd(), 'log'))
        outputs.append(os.path.join(os.getcwd(), 'plots'))
        outputs.append(os.path.join(os.getcwd(), 'simulation-state'))
        outputs.append(os.path.join(os.getcwd(), 'fission-rates'))
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


class TrackingTestHarness(TestHarness):
    """Specialized TestHarness for testing tracking."""
    
    def __init__(self):
        super(TrackingTestHarness, self).__init__()
        self.tracks = OrderedDict()
        self._result = ''

    def _segment_track(self, track, geometry):
        """Segments a given track over a given geometry and records the
           resulting segment information to a string"""

        # Segmentize a track in a geometry, recording the segments in a string
        geometry.segmentize(track)
        num_segments = track.getNumSegments()
        info = ' ' + str(num_segments) + '\n'
        for i in range(num_segments):
            info += str(i) + ': '
            segment = track.getSegment(i)
            info += str(round(segment._length, 8)) + ', '
            info += str(segment._region_id) + ', '
            info += str(segment._cmfd_surface_fwd) + ', '
            info += str(segment._cmfd_surface_bwd) + ', '
            info += str(segment._material.getName()) + ', '
            info += str(segment._material.getId()) + '\n'
        track.clearSegments()
        return info

    def _run_openmoc(self):
        """Segment tracks over the geometry and save the result to a string"""

        # Segmentize tracks over the geometry
        for track_name in self.tracks:
            self._result += track_name
            self._result += self._segment_track(self.tracks[track_name],
                                                self.input_set.geometry)

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""
        return self._result
