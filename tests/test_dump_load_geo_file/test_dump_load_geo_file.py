#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import NonUniformLatticeInput


class DumpAndLoadTestHarness(TestHarness):
    """An eigenvalue calculation with geometry load from and dump to a file."""

    def __init__(self, geo_type):
        super(DumpAndLoadTestHarness, self).__init__()
        self.input_set = NonUniformLatticeInput()
        self.num_polar = 4
        self.azim_spacing = self.spacing
        self.z_spacing = 0.5
        self.tolerance = 1E-4
        self.geo_type = geo_type


    def _create_geometry(self):
        """Initialize CMFD and add it to the Geometry."""

        self.input_set.geometry = openmoc.Geometry()
        if(self.geo_type ==1):
            # uniform lattice geometry file in old-versioned format
            self.input_set.geometry.loadFromFile('non-uniform-lattice1.geo', non_uniform_lattice=False)
            self.input_set.geometry.dumpToFile('log/dump1.geo', non_uniform_lattice=False)
        elif(self.geo_type ==2):
            # uniform lattice geometry file in new-versioned format
            self.input_set.geometry.loadFromFile('non-uniform-lattice2.geo', non_uniform_lattice=True)
            self.input_set.geometry.dumpToFile('log/dump2.geo', non_uniform_lattice=True)
        elif(self.geo_type ==3):
            # non-uniform lattice geometry file in new-versioned format
            self.input_set.geometry.loadFromFile('non-uniform-lattice3.geo', non_uniform_lattice=True)
            self.input_set.geometry.dumpToFile('log/dump3.geo', non_uniform_lattice=True)
        else:
            pass
        
        # Initialize CMFD
        cmfd = openmoc.Cmfd()
        cmfd.setCMFDRelaxationFactor(0.7)
        #cmfd.setLatticeStructure(17,17)
        cmfd.setWidths([[0.05,1.26,1.26,0.05], [0.05,1.26,1.26,0.05], [1.,1.5]])
        cmfd.setGroupStructure([[1,2,3],[4,5,6,7]])
        cmfd.setCentroidUpdateOn(False)

        # Add CMFD to the Geometry
        self.input_set.geometry.setCmfd(cmfd)
        
        

    def _create_trackgenerator(self):
        """Instantiate a TrackGenerator."""
        geometry = self.input_set.geometry
        geometry.initializeFlatSourceRegions()
        
        quad = openmoc.EqualWeightPolarQuad()
        quad.setNumPolarAngles(self.num_polar)
        self.track_generator = \
            openmoc.TrackGenerator3D(geometry, self.num_azim, self.num_polar,
                                     self.azim_spacing, self.z_spacing)
        self.track_generator.setQuadrature(quad)
        self.track_generator.setSegmentFormation(openmoc.OTF_STACKS)
        self.track_generator.setSegmentationZones([0.0,1.0,2.5])


    def _generate_tracks(self):
        """Generate Tracks and segments."""
        # Need to use more than 1 thread, and lose reproducibility of FSR
        # numbering, in order to have temporary tracks and segments array of
        # the correct size for the multi-threaded solver.
        self.track_generator.setNumThreads(self.num_threads)
        self.track_generator.generateTracks()


    def _create_solver(self):
        """Instantiate a CPULSSolver."""
        self.solver = openmoc.CPULSSolver(self.track_generator)
        self.solver.setNumThreads(self.num_threads)
        self.solver.setConvergenceThreshold(self.tolerance)


    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=False):
        """Digest info in the solver"""
        return super(DumpAndLoadTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)
        
    def _execute_test(self):
        """Build geometry, ray trace, run calculation."""
        self._setup()
        self._run_openmoc()

        
        

if __name__ == '__main__':
    try:
        harness1 = DumpAndLoadTestHarness(1)
        harness1._execute_test()
        harness2 = DumpAndLoadTestHarness(2)
        harness2._execute_test()
        harness3 = DumpAndLoadTestHarness(3)
        harness3._execute_test()
        results  = harness1._get_results()
        results += harness2._get_results()
        results += harness3._get_results()
        harness1._write_results(results)
        
        (harness1._opts, harness1._args) = harness1.parser.parse_args()
        if harness1._opts.update:
            harness1._overwrite_results()
        else :
            harness1._compare_results()
    finally:
        harness1._cleanup()
                    
