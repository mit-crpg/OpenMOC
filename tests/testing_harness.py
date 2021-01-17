from __future__ import print_function

import filecmp
import hashlib
import os
import shutil
import sys
import glob
import pickle
from collections import OrderedDict
from optparse import OptionParser
try:
    from PIL import Image, ImageOps
except ImportError as error:
    print(error.__class__.__name__ + ": " + str(error))

sys.path.insert(0, 'openmoc')
import openmoc
import openmoc.plotter
import openmoc.process

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


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
        # Add -f option to parser, as Jupyter forcedly provide this option
        self.parser.add_option('-f', dest='jupyter', action='store',default=None)

        self._opts = None
        self._args = None


    def _create_geometry(self):
        """Initialize the materials and geometry in the InputSet."""
        self.input_set.create_materials()
        self.input_set.create_geometry()

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        self.track_generator = \
            openmoc.TrackGenerator(geometry, self.num_azim, self.spacing)

    def _create_solver(self):
        """Instantiate a CPUSolver."""
        self.solver = openmoc.CPUSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)
        self.solver.setSolverMode(self.calculation_mode)

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

        # If running the test suite with MPI, only rank 0 should check the result
        try:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.Get_rank()
        except:
            rank = 0

        try:
            self._setup()
            self._run_openmoc()
            results = self._get_results()
            if rank == 0:
                self._write_results(results)
                self._compare_results()
        finally:
            if rank == 0:
                self._cleanup()

    def _update_results(self):
        """Update the results_true using the current version of OpenMOC."""

        # If running the test suite with MPI, only rank 0 should write the result
        try:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.Get_rank()
        except:
            rank = 0

        try:
            self._setup()
            self._run_openmoc()
            results = self._get_results()
            if rank == 0:
                self._write_results(results)
                self._overwrite_results()
        finally:
            if rank == 0:
                self._cleanup()

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue or fixed source calculation."""

        if self.solution_type == 'eigenvalue':
            self.solver.computeEigenvalue(self.max_iters, res_type=self.res_type)
        elif self.solution_type == 'flux':
            self.solver.computeFlux(self.max_iters)
        elif self.solution_type == 'source':
            self.solver.computeSource(self.max_iters, res_type=self.res_type)
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

        # For comparison of files with different line endings
        compare = (open('results_test.dat', 'r').read() ==
                   open('results_true.dat', 'r').read())
        if not compare:
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree.'

    def _cleanup(self):
        """Delete track, plot, etc. directories and test files."""

        # Create list of directories and/or files to remove
        outputs = [os.path.join(os.getcwd(), 'tracks')]
        outputs.append(os.path.join(os.getcwd(), 'log'))
        outputs.append(os.path.join(os.getcwd(), 'plots'))
        outputs.append(os.path.join(os.getcwd(), 'simulation-states'))
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
        self.zcoord = 0.0

    def _segment_track(self, track, geometry):
        """Segments a given track over a given geometry and records the
           resulting segment information to a string"""

        # Segmentize a track in a geometry, recording the segments in a string
        geometry.segmentize2D(track, self.zcoord)
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


class PlottingTestHarness(TestHarness):
    """Specialized TestHarness for testing plotting."""

    def __init__(self):
        super(PlottingTestHarness, self).__init__()
        self.figures = []

        # Use standardized default matplotlib rcparams
        rcparams = pickle.load(open('../rcparams.pkl', 'rb'))
        openmoc.plotter.matplotlib_rcparams = rcparams

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):

        # Store each each Matplotlib figure / PIL Image
        for i, fig in enumerate(self.figures):
            test_filename = 'test-{0}.png'.format(i)

            # Save the figure to a file
            if isinstance(fig, matplotlib.figure.Figure):
                fig.set_size_inches(4., 4.)
                fig.savefig(test_filename, bbox_inches='tight', dpi=100)
                plt.close(fig)
            else:
                fig.save(test_filename)

        return ''

    def _write_results(self, results_string):
        """Do nothing since the plots are created in _run_openmoc() method."""
        return

    def _overwrite_results(self):
        """Overwrite the reference images with the test images."""

        # Find all plot files
        outputs = glob.glob(os.path.join(os.getcwd(), 'test-*.png'))

        # Copy each test plot as a new reference plot
        for i in range(len(outputs)):
            shutil.copyfile('test-{0}.png'.format(i), 'true-{0}.png'.format(i))

    def _compare_results(self, max_distance=0.1):
        """Make sure the current results agree with the true standard."""

        # Loop over each Matplotlib figure / PIL Image and
        # compare to reference using Matplotlib fuzzy comparison
        for i, fig in enumerate(self.figures):

            # Open test image and resize to that of the true image with PIL
            img1 = Image.open('test-{0}.png'.format(i))
            img2 = Image.open('true-{0}.png'.format(i))
            img1 = ImageOps.fit(img1, img2.size, Image.ANTIALIAS)

            # Compute distance between each image in RGB space
            distance = self._compare_images(img1, img2)
            assert distance < max_distance, 'Results do not agree.'

    def _cleanup(self):
        """Delete plot PNG files."""

        # Find all test plot files
        outputs = glob.glob(os.path.join(os.getcwd(), 'test-*.png'))

        # Remove each plot file if it exists
        for i in range(len(outputs)):
            output = 'test-{0}.png'.format(i)
            if os.path.isfile(output):
                os.remove(output)
            elif os.path.isdir(output):
                shutil.rmtree(output)

        super(PlottingTestHarness, self)._cleanup()

    def _compare_images(self, img1, img2):
        """Compare two PIL images using in RGB space with pixel histograms."""

        # Extract RGBA data from each PIL Image
        rgba1 = np.array(img1)
        rgba2 = np.array(img2)

        # Compute histograms of each images pixels
        hr1, bins1 = np.histogram(rgba1[...,0], bins=256, density=True)
        hg1, bins1 = np.histogram(rgba1[...,1], bins=256, density=True)
        hb1, bins1 = np.histogram(rgba1[...,2], bins=256, density=True)
        hr2, bins2 = np.histogram(rgba2[...,0], bins=256, density=True)
        hg2, bins2 = np.histogram(rgba2[...,1], bins=256, density=True)
        hb2, bins2 = np.histogram(rgba2[...,2], bins=256, density=True)
        hist1 = np.array([hr1, hg1, hb1]).ravel()
        hist2 = np.array([hr2, hg2, hb2]).ravel()

        # Compute cartesian distance between histograms in RGB space
        diff = hist1 - hist2
        distance = np.sqrt(np.dot(diff, diff))
        return distance


class MultiSimTestHarness(TestHarness):
    """Specialized TestHarness for testing multi-simulation capabilities."""

    def __init__(self):
        super(MultiSimTestHarness, self).__init__()
        self.num_simulations = 3
        self.num_iters = []
        self.keffs = []

    def _run_openmoc(self):
        """Run multiple OpenMOC eigenvalue calculations."""

        for i in range(self.num_simulations):
            super(MultiSimTestHarness, self)._run_openmoc()
            self.num_iters.append(self.solver.getNumIterations())
            self.keffs.append(self.solver.getKeff())

    def _get_results(self, num_iterations=True, keff=True, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Return eigenvalues from each simulation into a string."""

        # Write out the iteration count and eigenvalues from each simulation
        outstr = ''
        for num_iters, keff in zip(self.num_iters, self.keffs):
            outstr += 'Iters: {0}\tkeff: {1:12.5E}\n'.format(num_iters, keff)

        return outstr
