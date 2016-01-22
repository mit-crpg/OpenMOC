#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import shutil
import re
import glob
import socket
from subprocess import call
from collections import OrderedDict
from optparse import OptionParser

# Command line parsing
parser = OptionParser()
parser.add_option('-j', '--parallel', dest='n_procs', default='1',
                  help="Number of parallel jobs.")
parser.add_option('-R', '--tests-regex', dest='regex_tests',
                  help="Run tests matching regular expression. \
                  Test names are the directories present in tests folder.\
                  This uses standard regex syntax to select tests.")
(options, args) = parser.parse_args()

# Define test data structure
tests = OrderedDict()

class Test(object):
    def __init__(self, name, debug=False):
        self.name = name
        self.debug = debug
        self.success = True
        self.msg = None
        self.skipped = False

        # FIXME: Set this up for different compilers???

    # Sets up build options for various tests. It is used both
    # in script and non-script modes
    def get_build_opts(self):
        build_str = ""
        if self.debug:
            build_str += "-Ddebug=ON "
        if self.coverage:
            build_str += "-Dcoverage=ON "
        self.build_opts = build_str
        return self.build_opts

    # FIXME: Run setup.py install
    def run_setup_install(self):
        """
        rc = call(self.cmake)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on setup.py'
        """

    # Runs make when in non-script mode
    def run_make(self):
        if not self.success:
            return

        # Default distutils build string
        make_list = ['python','setup.py', 'build_ext']

        # Use distutils to build OpenMOC
        rc = call(make_list)
        if rc != 0:
            self.success = False
            self.msg = 'Failed on setup.py'

    # Checks to see if file exists in PWD or PATH
    def check_compiler(self):
        result = False
        if os.path.isfile(self.fc):
            result = True
        for path in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(path, self.fc)):
                result = True
        if not result:
            self.msg = 'Compiler not found: {0}'.\
                       format((os.path.join(path, self.fc)))
            self.success = False

# Simple function to add a test to the global tests dictionary
def add_test(name, debug=False):
    tests.update({name: Test(name, debug)})

# List of all tests that may be run. User can add -C to command line to specify
# a subset of these configurations
add_test('normal')
add_test('debug', debug=True)

# Check to see if we should just print build configuration information to user
if options.list_build_configs:
    for key in tests:
        print('Configuration Name: {0}'.format(key))
        print('  Debug Flags:..........{0}'.format(tests[key].debug))
        print('  Coverage Test:........{0}\n'.format(tests[key].coverage))
    exit()

# Delete items of dictionary that don't match regular expression
if options.build_config is not None:
    for key in tests:
        if not re.search(options.build_config, key):
            del tests[key]

# Check if tests empty
if len(list(tests.keys())) == 0:
    print('No tests to run.')
    exit()

# Begin testing
shutil.rmtree('build', ignore_errors=True)
call(['./cleanup']) # removes all binary and hdf5 output files from tests
for key in iter(tests):
    test = tests[key]

    print('-'*(len(key) + 6))
    print(key + ' tests')
    print('-'*(len(key) + 6))
    sys.stdout.flush()

    # Verify fortran compiler exists
    test.check_compiler()
    if not test.success:
        continue

    if not test.success:
        continue

    # Get coverage command
    if test.coverage:
        test.find_coverage()
    if not test.success:
        continue

    # Run setup.py to configure build
    test.run_setup_install()

    # Go into build directory
    os.chdir('build')

    # Build OpenMC
    test.run_make()

    # Leave build directory
    os.chdir('..')

    # Copy over log file
    if script_mode:
        logfile = glob.glob('build/Testing/Temporary/LastTest_*.log')
    else:
        logfile = glob.glob('build/Testing/Temporary/LastTest.log')
    if len(logfile) > 0:
        logfilename = os.path.split(logfile[0])[1]
        logfilename = os.path.splitext(logfilename)[0]
        logfilename = logfilename + '_{0}.log'.format(test.name)
        shutil.copy(logfile[0], logfilename)

    # Clear build directory and remove binary and hdf5 files
    shutil.rmtree('build', ignore_errors=True)
    if script_mode:
        os.remove('ctestscript.run')
    call(['./cleanup'])

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
