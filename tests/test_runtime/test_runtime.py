#!/usr/bin/env python

import os
import sys
import numpy as np
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
from openmoc.log import py_printf

class RuntimeTestHarness(TestHarness):
    '''An eigenvalue calculation in a 3D lattice with 7-group C5G7 data.'''

    def __init__(self):
        super(RuntimeTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(3)

    def _create_trackgenerator(self):
        '''Dump the geometry, to be loaded through runtime class'''

        # Geometry should be dumped before FSRs are initialized
        # Dump the geometry
        self.input_set.geometry.dumpToFile('geometry_file.geo')

        # Get rid of the geometry
        self.input_set.geometry = None

    def _create_solver(self):
        pass

    def _generate_tracks(self):
        pass

    def _run_openmoc(self):

        runtime = openmoc.RuntimeParameters()

        # Display help message
        runtime.setRuntimeParameters(["--help".encode('utf8')])
        
        string_input = ['-debug', '1',
        '-log_level', 'NORMAL',
        '-domain_decompose', '2,2,2',
        '-num_domain_modules', '1,1,1',
        '-num_threads', '1',
        '-log_filename', 'test_problem.log',
        '-geo_filename', 'geometry_file.geo',
        '-azim_spacing', '0.22 ',
        '-num_azim', '4',
        '-polar_spacing', '0.8',
        '-num_polar', '6',
        '-seg_zones', '-5.0,5.0',
        '-segmentation_type', '3',
        '-quadraturetype', '2',
        '-CMFD_group_structure', '1/2,3,4/5,6,7',
        '-CMFD_lattice', '2,3,3',
        '-widths_x', '1,2*1,1',
        '-widths_y', '1,2*1,1',
        '-widths_z', '3.5,2.5*2,1.5',
        '-CMFD_flux_update_on', '1',
        '-knearest', '2',
        '-CMFD_centroid_update_on', '1',
        '-use_axial_interpolation', '2',
        '-SOR_factor', '1.5',
        '-CMFD_relaxation_factor', '0.7',
        '-ls_solver', '1',
        '-max_iters', '15',
        '-MOC_src_residual_type', '1',
        '-MOC_src_tolerance', '1.0E-2',
        '-output_mesh_lattice', '-output_mesh_lattice', ' 5,5,9', ' -output_type', ' 0',
        '-output_mesh_lattice', '-output_mesh_lattice', ' 5,5,9', ' -output_type', ' 1',
        '-non_uniform_output', '1.26*3/1*3/4.*3/-1.,1.,-1.', ' -output_type 1  ',
        '-verbose_report', '1',
        '-time_report', '1']
        string_input = [s.encode('utf8') for s in string_input]

        runtime.setRuntimeParameters(string_input)

        # Define simulation parameters
        num_threads = runtime._num_threads;

        # Set logging information
        if (runtime._log_filename):
            openmoc.set_log_filename(runtime._log_filename);
        openmoc.set_log_level(runtime._log_level);
        openmoc.set_line_length(120);

        py_printf('NORMAL', 'Geometry file = %s', runtime._geo_filename)
        py_printf('NORMAL', 'Azimuthal spacing = %f', runtime._azim_spacing)
        py_printf('NORMAL', 'Azimuthal angles = %d', runtime._num_azim)
        py_printf('NORMAL', 'Polar spacing = %f', runtime._polar_spacing)
        py_printf('NORMAL', 'Polar angles = %d', runtime._num_polar)

        # Create the geometry
        py_printf('NORMAL', 'Creating geometry...');
        geometry = openmoc.Geometry()
        self.input_set.geometry = geometry
        if (not runtime._geo_filename):
            py_printf('ERROR', 'No geometry file is provided');
        geometry.loadFromFile(runtime._geo_filename)
        if False: #FIXME
            geometry.setDomainDecomposition(runtime._NDx, runtime._NDy, runtime._NDz, 
                                            MPI_COMM_WORLD)

        geometry.setNumDomainModules(runtime._NMx, runtime._NMy, runtime._NMz)

        if ((runtime._NCx >0 and runtime._NCy > 0 and runtime._NCz > 0) or
            (not runtime._cell_widths_x.empty() and not runtime._cell_widths_y.empty() and
             not runtime._cell_widths_z.empty())):

            # Create CMFD mesh
            py_printf('NORMAL', 'Creating CMFD mesh...')
            cmfd = openmoc.Cmfd()
            cmfd.setSORRelaxationFactor(runtime._SOR_factor)
            cmfd.setCMFDRelaxationFactor(runtime._CMFD_relaxation_factor)
            if (runtime._cell_widths_x.empty() or runtime._cell_widths_y.empty() or
                runtime._cell_widths_z.empty()):
                cmfd.setLatticeStructure(runtime._NCx, runtime._NCy, runtime._NCz)
            else:
                cmfd_widths = [runtime._cell_widths_x, runtime._cell_widths_y, runtime._cell_widths_z]
                cmfd.setWidths(cmfd_widths);

            if (not runtime._CMFD_group_structure.empty()):
                cmfd.setGroupStructure(runtime._CMFD_group_structure)
            cmfd.setKNearest(runtime._knearest)
            cmfd.setCentroidUpdateOn(runtime._CMFD_centroid_update_on)
            cmfd.useAxialInterpolation(runtime._use_axial_interpolation)

            geometry.setCmfd(cmfd)

        geometry.initializeFlatSourceRegions()

        # Initialize track generator and generate tracks
        py_printf('NORMAL', 'Initializing the track generator...');
        if runtime._quadraturetype == 0:
            quad = openmoc.TYPolarQuad()
        if runtime._quadraturetype == 1:
            quad = openmoc.LeonardPolarQuad()
        if runtime._quadraturetype == 2:
            quad = openmoc.GLPolarQuad()
        if runtime._quadraturetype == 3:
            quad = openmoc.EqualWeightPolarQuad()
        if runtime._quadraturetype == 4:
            quad = openmoc.EqualAnglePolarQuad()

        quad.setNumAzimAngles(runtime._num_azim);
        quad.setNumPolarAngles(runtime._num_polar);
        track_generator = openmoc.TrackGenerator3D(geometry, runtime._num_azim, 
                                   runtime._num_polar, runtime._azim_spacing,
                                   runtime._polar_spacing)
        track_generator.setNumThreads(num_threads)
        track_generator.setQuadrature(quad)
        track_generator.setSegmentFormation(runtime._segmentation_type)
        if(len(runtime._seg_zones) > 0):
            track_generator.setSegmentationZones(runtime._seg_zones)
        track_generator.generateTracks()
        self.track_generator = track_generator

        # Initialize solver and run simulation
        if (runtime._linear_solver):
            self.solver = openmoc.CPULSSolver(track_generator)
        else:
            self.solver = openmoc.CPUSolver(track_generator)
        if (runtime._verbose_report):
            self.solver.setVerboseIterationReport()
        self.solver.setNumThreads(num_threads)
        self.solver.setConvergenceThreshold(runtime._tolerance)
        self.solver.computeEigenvalue(runtime._max_iters, runtime._MOC_src_residual_type)

        if (runtime._time_report):
            self.solver.printTimerReport()

        # Extract reaction rates
        my_rank = 0
        #if False: #FIXME
            #MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        rxtype = {'FISSION_RX', 'TOTAL_RX', 'ABSORPTION_RX', 'FLUX_RX'}

        # Consider using the mesh output as the test result



    def _get_results(self, num_iters=True, keff=True, fluxes=False,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=False):
        '''Digest info in the solver'''
        return super(RuntimeTestHarness, self)._get_results(
                num_iters=num_iters, keff=keff, fluxes=fluxes,
                num_fsrs=num_fsrs, num_tracks=num_tracks,
                num_segments=num_segments, hash_output=hash_output)

if __name__ == '__main__':
    harness = RuntimeTestHarness()
    harness.main()
