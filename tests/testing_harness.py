import filecmp
import hashlib
import os
import shutil

# from input_set import InputSet


class TestHarness(object):
    """General class for running OpenMOC regression tests."""
    def __init__(self, statepoint_name, tallies_present=False):
        # FIXME: Use OpenMOC Options module here???
        self._input_set = InputSet()
        self._opts = None
        # self._solver = None
        self.parser.add_option('--update', dest='update', action='store_true',
                               default=False)

    def main(self):
        """Accept command line arguments and either run or update tests."""
        self._parse_args()
        if self._opts.update:
            self.update_results()
        else:
            self.execute_test()

    def execute_test(self):
        """Run OpenMOC with the appropriate arguments and check the outputs."""
        try:
            self._build_inputs()
            self._run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update the results_true using the current version of OpenMOC."""
        try:
            self._build_inputs()
            self._run_openmoc()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()

    def _build_inputs(self):
        """Build OpenMOC TrackGenerator and Solver from InputSet."""
        self._input_set.build_default_materials_and_geometry()
        self._input_set.build_default_settings()
        self._input_set.export()

    def _get_results(self, hash_output=False):
        """Digest info from solver and return as a string."""

        # FIXME: Write out keff
        # FIXME: Write out FSR fiss rates, fluxes? Is this reproducible?

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
        """Delete track, simulation state and plot files."""
        shutil.rmtree('/tracks')
        shutil.rmtree('/simulation-state')
        shutil.rmtree('/plots')


class HashedTestHarness(TestHarness):
    """Specialized TestHarness that hashes the results."""
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        return TestHarness._get_results(self, True)


class PlotTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMOC plotting tests."""
    def __init__(self, plot_names):
        self._plot_names = plot_names
        self._opts = None
        self._args = None

    def _run_openmoc(self):
        assert returncode == 0, 'OpenMOC did not exit successfully.'

    def _cleanup(self):
        TestHarness._cleanup(self)
        output = glob.glob(os.path.join(os.getcwd(), '*.ppm'))
        for f in output:
            if os.path.exists(f):
                os.remove(f)

    def _get_results(self):
        """Return a string hash of the plot files."""
        # Find the plot files
        # FIXME: Does this work???
        plot_files = glob.glob(os.path.join(os.getcwd() + '/plots', '*.png'))
        print(plot_files)

        # Read the plot files
        outstr = bytes()
        for fname in sorted(plot_files):
            with open(fname, 'rb') as fh:
                outstr += fh.read()

        # Hash the information and return
        sha512 = hashlib.sha512()
        sha512.update(outstr)
        outstr = sha512.hexdigest()

        return outstr
