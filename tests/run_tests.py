#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import shutil
import re
import glob
import copy
import subprocess
from collections import OrderedDict
from optparse import OptionParser
from distutils.spawn import find_executable

# Command line parsing
parser = OptionParser()
parser.add_option('-j', '--parallel', dest='n_procs', default='1',
                  help="Number of parallel jobs.")
parser.add_option('-v', '--verbose', action="store_true", dest='verbose',
                  default=False, help="Make CTest verbose.")
parser.add_option('-R', '--tests-regex', dest='regex_tests',
                  help="Run tests matching regular expression. "
                  "Test names are the directories present in tests folder. "
                  "This uses standard regex syntax to select tests.")
parser.add_option('-C', '--build-config', dest='build_config',
                  help="Build configurations matching regular expression. "
                        "Specific build configurations can be printed out with "
                        "optional argument -p, --print. This uses standard "
                        "regex syntax to select build configurations.")
parser.add_option('-c', '--coverage', action="store_true", dest='coverage',
                   default=False, help="Run tests with coverage.py to output "
                  "Python code coverage.")
parser.add_option('-l', '--list', action="store_true",
                  dest="list_build_configs", default=False,
                  help="List out build configurations.")
parser.add_option('-m', '--with-mpi', action="store_true",
                  dest="with_mpi", default=False,
                  help="Run MPI build configuration")
(options, args) = parser.parse_args()

# Default build options
FP = 'double'

# Define test data structure
tests = OrderedDict()

class Test(object):

    # A class attribute to cache the setup install commands from
    # previous Tests, if any (this helps eliminate redundant builds)
    _setup_cmd = []

    def __init__(self, name, cc='gcc', num_threads=1, debug=False, 
                 coverage=False):
        self.name = name
        self.cc = cc
        self.fp = FP
        self.num_threads = num_threads
        self.debug = debug
        self.coverage = coverage
        self.success = True
        self.msg = None

    def run_setup_install(self):
        """Install OpenMOC with distutils"""

        setup_cmd = [sys.executable, 'setup.py', 'install']
        setup_cmd += ['--install-purelib=tests/openmoc']
        setup_cmd += ['--cc={0}'.format(self.cc), '--fp={0}'.format(self.fp)]
        if self.debug:
            setup_cmd += ['--debug-mode']
        if self.coverage:
            setup_cmd += ['--coverage-mode']

        # Run setup.py if it was not run for the previous Test
        if setup_cmd != Test._setup_cmd:
            rc = subprocess.call(setup_cmd)
            rc = subprocess.call(setup_cmd)

            # Check for error code
            if rc != 0:
                self.success = False
                self.msg = 'Failed on setup.py'
            # Cache the setup install command for the next Test
            else:
                Test._setup_cmd = setup_cmd

    def run_cmake(self):
        """Run CMake to create CTest script"""

        cmake_cmd = ['cmake', '-H..', '-Bbuild']

        # Run tests with coverage option or not
        if options.coverage:
            # Python API coverage testing
            coverage_exe = find_executable('coverage')
            if not coverage_exe:
                raise RuntimeError('Coverage executable was not found.')

            cmake_cmd += ['-DPYTHON_EXECUTABLE=' + coverage_exe]
            cmake_cmd += ['-DCOVERAGE=True']
            cmake_cmd += ['-DOMIT=test*']
        else:
            cmake_cmd += ['-DPYTHON_EXECUTABLE=' + sys.executable]

        # Run CMake
        rc = subprocess.call(cmake_cmd)

        # Check for error code
        if rc != 0:
            self.success = False
            self.msg = 'Failed on cmake.'

    def run_ctests(self):
        """Run CTest on all tests"""
        if not self.success:
            return

        os.environ['OMP_NUM_THREADS'] = str(self.num_threads)

        # Default CTest string
        ctest_cmd = ['ctest']

        # Check for verbose
        if options.verbose:
            ctest_cmd += ['--verbose']

        # Check for running tests in parallel
        if options.n_procs:
            ctest_cmd += ['-j', str(options.n_procs)]

        # Check for subset of tests
        if options.regex_tests:
            ctest_cmd += ['-R', str(options.regex_tests)]

        # Avoid MPI tests when running a non-mpi build
        if 'mpi' not in self.name:
            ctest_cmd += ['-E', 'mpi']
        #FIXME Skip runtime test that stopped working when adding --help
        else:
            ctest_cmd += ['-E', 'runtime']

        # Run CTest
        rc = subprocess.call(ctest_cmd)

        # Check for error code
        if rc != 0:
            self.success = False
            self.msg = 'Failed on testing.'


# Simple function to add a test to the global tests dictionary
def add_test(name, cc='gcc', num_threads=1, debug=False, ):
    tests.update({name: Test(name, cc, num_threads, debug, options.coverage)})

# Look at environment variables to see which tests should be run
omp = (os.environ.get('OMP') == 'y')
mpi = (os.environ.get('MPI') == 'y')

# List of all tests that may be run. User can add -C to command line to specify
# a subset of these configurations
if not omp and not mpi:
    add_test('normal-gcc', cc='gcc', num_threads=1)
elif omp and not mpi:
    add_test('normal-openmp-gcc', cc='gcc', num_threads=4)
elif (omp and mpi) or options.with_mpi:
    add_test('mpi-openmp-gcc', cc='mpicc', num_threads=2)
#add_test('normal-gcc-debug', cc='gcc', num_threads=1, debug=True)
#add_test('normal-openmp-gcc-debug', cc='gcc', num_threads=4, debug=True)
#add_test('normal-icpc', cc='icpc', num_threads=1)
#add_test('normal-openmp-icpc', cc='icpc', num_threads=4)
#add_test('normal-clang', cc='clang', num_threads=1)
#add_test('normal-openmp-clang', cc='clang', num_threads=4)

# Check to see if we should just print build configuration information to user
if options.list_build_configs:
    for key in tests:
        print('Configuration Name: {0}'.format(key))
    exit()

# Delete items of dictionary that don't match regular expression
if options.build_config is not None:
    tests_copy = copy.deepcopy(tests)
    for key in tests:
        if not re.search(options.build_config, key):
            del tests_copy[key]
    tests = tests_copy

# Check if tests empty
if len(list(tests.keys())) == 0:
    print('No tests to run.')
    exit()

# Removes all binary track and output files from tests
shutil.rmtree('build', ignore_errors=True)
subprocess.call(['./cleanup'])

# Run each Test in sequence
for key in iter(tests):
    test = tests[key]

    # Print header for this test
    print('-'*(len(key) + 6))
    print(key + ' tests')
    print('-'*(len(key) + 6))
    sys.stdout.flush()

    # Run CMake to setup CTest
    test.run_cmake()

    # Go into main OpenMOC directory
    os.chdir('..')

    # Run setup.py to build and install OpenMOC
    test.run_setup_install()

    # Go into build directory
    os.chdir('tests/build')

    # Run CTest
    test.run_ctests()

    # Leave build directory
    os.chdir('..')

    # Copy over log file
    logfile = glob.glob('build/Testing/Temporary/LastTest.log')
    if len(logfile) > 0:
        logfilename = os.path.split(logfile[0])[1]
        logfilename = os.path.splitext(logfilename)[0]
        logfilename = logfilename + '_{0}.log'.format(test.name)
        shutil.copy(logfile[0], logfilename)

# Combine all coverage files
if options.coverage:
    os.system('coverage combine test*/.coverage unit_tests/.coverage')

# Clear build directory and remove binary and hdf5 files
shutil.rmtree('build', ignore_errors=True)
if not options.coverage:
    shutil.rmtree('openmoc', ignore_errors=True)
subprocess.call(['./cleanup'])

# Print out summary of results
print('\n' + '='*54)
print('Summary of Compilation Option Testing:\n')

if sys.stdout.isatty():
    OK = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
else:
    OK = ''
    FAIL = ''
    ENDC = ''
    BOLD = ''

return_code = 0

for test in tests:
    print(test + '.'*(50 - len(test)), end='')
    if tests[test].success:
        print(BOLD + OK + '[OK]' + ENDC)
    else:
        print(BOLD + FAIL + '[FAILED]' + ENDC)
        print(' '*len(test)+tests[test].msg)
        return_code = 1

sys.exit(return_code)
