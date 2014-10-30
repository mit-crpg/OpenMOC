## Works when run from commandline using python test-options.py

## DOES NOT work from IDLE (Python crashes when importing options)
## DOES NOT work with nosetests from command line (module does not have attribute 'test-options')


import unittest
#import nose
import sys
from sys import path
import openmoc.options as options
import multiprocessing


class TestDefaultInit(unittest.TestCase):
    
    ## Tests the default values set when creating an Options instance (when nothing is specified in the command line)

    @classmethod
    def setUpClass(cls):
        cls._default_options = options.Options()

    def test_default_num_azim(self):
        self.assertEqual(self._default_options.getNumAzimAngles(), 4)

    def test_default_track_spacing(self):
        self.assertEqual(self._default_options.getTrackSpacing(), 0.1)

    def test_default_max_iters(self):
        self.assertEqual(self._default_options.getMaxIterations(), 1000)

    def test_default_tolerance(self):
        self.assertEqual(self._default_options.getTolerance(), .00001)

    def test_default_num_omp_threads(self):
        self.assertEqual(self._default_options.getNumThreads(), multiprocessing.cpu_count())

    def test_default_num_thread_blocks(self):
        self.assertEqual(self._default_options.getNumThreadBlocks(), 64)

    def test_default_num_gpu_threads(self):
        self.assertEqual(self._default_options.getNumThreadsPerBlock(), 64)

    def test_default_cmfd_acceleration(self):
        self.assertFalse(self._default_options.getCmfdAcceleration())

    def test_default_cmfd_relax_factor(self):
        self.assertEqual(self._default_options.getCmfdRelaxationFactor(), 0.6)

    def test_default_cmfd_mesh_level(self):
        self.assertEqual(self._default_options.getCmfdMeshLevel(), -1)

    @classmethod
    def tearDownClass(cls):
        cls._default_options = None

class TestCustomInit(unittest.TestCase):

    # Uses simulated command line input, test parseArguments
    
    @classmethod
    def setUpClass(cls):

        ## directly changing sys.argv to reflect the following command line:

        ## --num-azim 5 --track-spacing 0.4 --max-iters 2000 --tolerance 0.00005 --num-omp-threads 2 --num-thread-blocks 72
        ## --num-gpu-threads 80 --relax-factor 0.7 --mesh-level -2

        sys.argv = ['test-options.py']
        sys.argv.extend(["--num-azim", "5", '--track-spacing', '0.4', '--max-iters', '2000', '--tolerance', '0.00005', '--num-omp-threads', '2',
                         '--num-thread-blocks', '72', '--num-gpu-threads', '80', '--relax-factor', '0.7', '--mesh-level', '-2'])
        
        cls._options = options.Options()

        ## parseArguments is part of initializing Options(), so it will get the values from simulated command line

        
    def test_command_num_azim(self):
        self.assertEqual(self._options.getNumAzimAngles(), 5)

    def test_command_track_spacing(self):
        self.assertEqual(self._options.getTrackSpacing(), 0.4)

    def test_command_max_iters(self):
        self.assertEqual(self._options.getMaxIterations(), 2000)

    def test_command_tolerance(self):
        self.assertEqual(self._options.getTolerance(), .00005)

    def test_command_num_omp_threads(self):
        self.assertEqual(self._options.getNumThreads(), 2)

    def test_command_num_thread_blocks(self):
        self.assertEqual(self._options.getNumThreadBlocks(), 72)

    def test_command_num_gpu_threads(self):
        self.assertEqual(self._options.getNumThreadsPerBlock(), 80)

    def test_command_cmfd_acceleration(self):
        self.assertFalse(self._options.getCmfdAcceleration())

    def test_command_cmfd_relax_factor(self):
        self.assertEqual(self._options.getCmfdRelaxationFactor(), 0.7)

    def test_command_cmfd_mesh_level(self):
        self.assertEqual(self._options.getCmfdMeshLevel(), -2)

    @classmethod
    def tearDownClass(cls):
        cls._options = None
        sys.argv = ['test-options.py']
    

TestOptionsSuite = unittest.TestLoader().loadTestsFromTestCase(TestDefaultInit)
TestOptionsSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCustomInit))

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(TestOptionsSuite)
