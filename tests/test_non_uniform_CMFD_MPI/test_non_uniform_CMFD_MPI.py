#!/usr/bin/env python

import os
import sys
import glob
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
import openmoc
from openmoc.log import py_printf
import subprocess


class TestNonUniformCMFDMPITestHarness(TestHarness):
    """A test of non-uniform CMFD OpenMOC running in parallel."""

    def __init__(self):
        super(TestNonUniformCMFDMPITestHarness, self).__init__()

    def _create_geometry(self):
        pass

    def _create_trackgenerator(self):
        pass

    def _generate_tracks(self):
        pass
        
    def _creat_solver(self):
        pass

    def _run_openmoc(self):
      pwd = os.getcwd()
      cmd = 'cp ../../profile/models/run_time_standard/run_time_standard.cpp ' \
            + pwd
      cmd += ' && ' + 'cp ../../profile/models/run_time_standard/' + \
            'non-uniform-lattice.geo ' + pwd
      rc = subprocess.call(cmd, shell=True)
      
       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to copy run_time_standard.cpp or \
              non-uniform-lattice.geo.\n')
        os._exit(0)
      
      
      os.chdir("../../profile/")
      case = pwd + '/run_time_standard.cpp'
      cmd = 'make clean case=' + case + ' && make -j 10 OPTIMIZE=yes \
            DEBUG=no case=' + case
      rc = subprocess.call(cmd, shell=True)
      
       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to build the execute file. Check the Makefile \
              or source code\n')
        os._exit(0)

      
      
      os.chdir(pwd)
      
      # 1. Running in uniform CMFD with setLatticeStructure, without MPI
      cmd = """mpiexec -n 1 ./run_time_standard -domain_decompose 1,1,1 \
      -CMFD_lattice 2,2,2 -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #1. \
              Check the parameters\n')
        os._exit(0)
      
      # 2. Running in uniform CMFD with setLatticeStructure, with MPI
      cmd = """mpiexec -n 8 ./run_time_standard -domain_decompose 2,2,2 \
      -CMFD_lattice 2,2,2 -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #2. \
              Check the parameters\n')
        os._exit(0)
      
      # 3. Running in uniform CMFD with set uniform  XYZ widths, without MPI
      cmd = """mpiexec -n 1 ./run_time_standard -domain_decompose 1,1,1 \
      -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -widths_x 1.31*2 -widths_y 1.31*2 -widths_z 1.25*2 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #3. \
              Check the parameters\n')
        os._exit(0)
      
      # 4. Running in uniform CMFD with set uniform  XYZ widths, with MPI
      cmd = """mpiexec -n 8 ./run_time_standard -domain_decompose 2,2,2 \
      -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -widths_x 1.31*2 -widths_y 1.31*2 -widths_z 1.25*2 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #4. \
              Check the parameters\n')
        os._exit(0)
      
      # 5. Running in non-uniform CMFD, without MPI
      cmd = """mpiexec -n 1 ./run_time_standard -domain_decompose 1,1,1 \
      -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -widths_x 0.05,1.26*2,0.05 -widths_y 0.05,1.26*2,0.05 -widths_z 1.25*2 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #5. \
              Check the parameters\n')
        os._exit(0)
      
      # 6. Running in non-uniform CMFD, with MPI
      cmd = """mpiexec -n 8 ./run_time_standard -domain_decompose 2,2,2 \
      -SOR_factor 1.5 -CMFD_relaxation_factor 0.7 \
      -CMFD_flux_update_on 1 \
      -widths_x 0.05,1.26*2,0.05 -widths_y 0.05,1.26*2,0.05 -widths_z 1.25*2 \
      -seg_zones 0.,2.5  -max_iters 1000 \
      -log_filename non-uniform-lattice_test.log \
      -geo_file_name non-uniform-lattice.geo"""
      rc = subprocess.call(cmd, shell=True)                

       # Check for error code
      if rc != 0:
        print('\nERROR: Failed to run the test problem #6. \
              Check the parameters\n')
        os._exit(0)
      
    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
      """Digest info in the log file and return as a string."""

      # Find the log filename with the time and date
      logfilename = glob.glob('log/non-uniform-lattice_test.log')
      
      with open(logfilename[0], 'r') as f:
        lines = f.readlines()
          
      with open(logfilename[0], 'w') as f_w:
        for line in lines:
          if "Current local time and date" in line:
            continue          
          if " sec" in line:
            continue
          i = f_w.write(line)
              
      # Read the file into a list of strings for each line
      with open(logfilename[0], 'r') as myfile:
        lines = myfile.readlines()
          
      cmd = 'rm -rf ./log run_time_standard* non-uniform-lattice.geo'
      rc = subprocess.call(cmd, shell=True) 
      
      # Concatenate all strings in the file into a single string
      outstr = ''.join(lines)
      return outstr


if __name__ == '__main__':
    harness = TestNonUniformCMFDMPITestHarness()
    harness.main()
