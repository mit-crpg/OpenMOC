#!/usr/bin/env python

import os
import sys
import glob
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
import openmoc
from openmoc.log import py_printf


class LoggerTestHarness(TestHarness):
    """A test of OpenMOC's logger module."""

    def __init__(self):
        super(LoggerTestHarness, self).__init__()

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass

    def _run_openmoc(self):
        """Print a variety of log messages to a log file."""

        # Set a log level which precludes some messages from being printed
        openmoc.set_log_level('NORMAL')

        # Print messages using the pure C implementation
        openmoc.log_printf(openmoc.DEBUG, 'This is a debug message')
        openmoc.log_printf(openmoc.INFO, 'This is an info message')
        openmoc.log_printf(openmoc.NORMAL, 'This is a normal message')
        openmoc.log_printf(openmoc.SEPARATOR, 'This is a separator message')
        openmoc.log_printf(openmoc.HEADER, 'This is a header message')
        openmoc.log_printf(openmoc.TITLE, 'This is a title message')
        openmoc.log_printf(openmoc.WARNING, 'This is a warning message')
        openmoc.log_printf(openmoc.CRITICAL, 'This is a critical message')
        openmoc.log_printf(openmoc.RESULT, 'This is a result message')

        # Print messages using the Python-wrapped version
        py_printf('DEBUG', 'This is a debug message')
        py_printf('INFO', 'This is an info message')
        py_printf('NORMAL', 'This is a normal message')
        py_printf('SEPARATOR', 'This is a separator message')
        py_printf('HEADER', 'This is a header message')
        py_printf('TITLE', 'This is a title message')
        py_printf('WARNING', 'This is a warning message')
        py_printf('CRITICAL', 'This is a critical message')
        py_printf('RESULT', 'This is a result message: %d', 5)

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
        """Digest info in the log file and return as a string."""

        # Find the log filename with the time and date
        logfilename = glob.glob('log/openmoc-*')

        # Read the file into a list of strings for each line
        with open(logfilename[0], 'r') as myfile:
            lines = myfile.readlines()
        
        # Concatenate all strings in the file into a single string
        # Exclude the first line which is the time and date
        outstr = ''.join(lines[1:])
        return outstr


if __name__ == '__main__':
    harness = LoggerTestHarness()
    harness.main()
