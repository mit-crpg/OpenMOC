#!/usr/bin/env python

import unittest
import os
import getopt
import sys
import multiprocessing
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc.options

# This unit test checks the openmoc.options Options class. It makes sure the
# default parameters stay the same, and more importantly checks that the parser
# keeps working.

class TestDefaultInit(unittest.TestCase):

    # creates a new Options() every test
    def setUp(self):
        self._default_options = openmoc.options.Options()

    def test_default_num_azim(self):
        self.assertEqual(self._default_options._num_azim, 4)

    def test_default_num_azim(self):
        self.assertEqual(self._default_options._num_polar, 6)

    def test_default_polar_spacing(self):
        self.assertEqual(self._default_options._polar_spacing, 1.5)

    def test_default_track_spacing(self):
        self.assertEqual(self._default_options._azim_spacing, 0.1)

    def test_default_max_iters(self):
        self.assertEqual(self._default_options._max_iters, 1000)

    def test_default_tolerance(self):
        self.assertEqual(self._default_options._tolerance, .00001)

    def test_default_num_omp_threads(self):
        try:
            default_num_threads = int(os.environ["OMP_NUM_THREADS"])
        except KeyError:
            default_num_threads = multiprocessing.cpu_count()
        self.assertEqual(self._default_options._num_omp_threads,
             default_num_threads)

    def test_default_num_thread_blocks(self):
        self.assertEqual(self._default_options._num_thread_blocks, 64)

    def test_default_num_gpu_threads(self):
        self.assertEqual(self._default_options._num_threads_per_block, 64)

    def tearDown(self):
        self._default_options = None

class TestCustomInit(unittest.TestCase):

    def setUp(self):
        sys.argv = ["dummy.py", "-a", "8", "-s", "0.2", "-p", "10", "-l",
                    "1.2", "-i", "100", "-c", "0.1", "-t", "4", "-b", "32",
                    "-g", "32"]
        self._custom_options = openmoc.options.Options()

    def test_custom_num_azim(self):
        self.assertEqual(self._custom_options._num_azim, 8)

    def test_custom_num_azim(self):
        self.assertEqual(self._custom_options._num_polar, 10)

    def test_custom_polar_spacing(self):
        self.assertEqual(self._custom_options._polar_spacing, 1.2)

    def test_custom_track_spacing(self):
        self.assertEqual(self._custom_options._azim_spacing, 0.2)

    def test_custom_max_iters(self):
        self.assertEqual(self._custom_options._max_iters, 100)

    def test_custom_tolerance(self):
        self.assertEqual(self._custom_options._tolerance, .1)

    def test_custom_num_omp_threads(self):
        self.assertEqual(self._custom_options._num_omp_threads, 4)

    def test_custom_num_thread_blocks(self):
        self.assertEqual(self._custom_options._num_thread_blocks, 32)

    def test_custom_num_gpu_threads(self):
        self.assertEqual(self._custom_options._num_threads_per_block, 32)

suite = unittest.TestLoader().loadTestsFromTestCase(TestDefaultInit)
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCustomInit))

unittest.TextTestRunner(verbosity=2).run(suite)
