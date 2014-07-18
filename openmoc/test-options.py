## Works when run from commandline using python test-options.py

## DOES NOT work from IDLE (Python crashes when importing options)
## DOES NOT work with nosetests from command line (module does not have attribute 'test-options')

import unittest
import nose
from options import *

class TestDefaultInit(unittest.TestCase):
    
    # creates a new Options() every test... can it just be done once?

    def setUp(self):
        self._default_options = Options()

    def test_default_num_azim(self):
        self.assertEqual(self._default_options._num_azim, 4)

    def test_default_track_spacing(self):
        self.assertEqual(self._default_options._track_spacing, 0.1)

    def test_default_max_iters(self):
        self.assertEqual(self._default_options._max_iters, 1000)

    def test_default_tolerance(self):
        self.assertEqual(self._default_options._tolerance, .00001)

    def test_default_num_omp_threads(self):
        self.assertEqual(self._default_options._num_omp_threads, 1)

    def test_default_num_thread_blocks(self):
        self.assertEqual(self._default_options._num_thread_blocks, 64)

    def test_default_num_gpu_threads(self):
        self.assertEqual(self._default_options._num_gpu_threads, 64)

    def test_default_cmfd_acceleration(self):
        self.assertFalse(self._default_options._use_cmfd_acceleration)

    def test_default_cmfd_relax_factor(self):
        self.assertEqual(self._default_options._cmfd_relax_factor, 0.6)

    def test_default_cmfd_mesh_level(self):
        self.assertEqual(self._default_options._cmfd_mesh_level, -1)

    def tearDown(self):
        self._default_options = None

class TestCustomInit(unittest.TestCase):

    def setUp(self):
        self._options = Options()
        opts, args = getopt.getopt('REPLACE THIS WITH THE ARGS',
                                 'hfa:s:i:c:t:b:g:r:l:',
                                 ['help',
                                  'num-azim=',
                                  'track-spacing=',
                                  'tolerance=',
                                  'max-iters=',
                                  'num-omp-threads=',
                                  'num-thread-blocks=',
                                  'num-gpu-threads=',
                                  'relax-factor=',
                                  'cmfd-acceleration',
                                  'mesh-level='])
    

suite = unittest.TestLoader().loadTestsFromTestCase(TestDefaultInit)
#suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCustomInit))

unittest.TextTestRunner(verbosity=2).run(suite)        
