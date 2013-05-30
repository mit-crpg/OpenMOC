###############################################################################
#                                 User Options
###############################################################################

# Name of the Python package to be created
package_name = 'openmoc'

# Supported C++ compilers: 'gcc', 'icpc', 'all'
cpp_compilers = ['all']

# Supported floating point precision: 'single', 'double', 'all'
fp_precision = ['all']

# Use CUDA set to True or False
with_cuda = True



###############################################################################
#                                 Source Code
###############################################################################

sources = {}

sources['c++'] = ['openmoc/openmoc.i',
                  'openmoc/src/host/Cell.cpp',
                  'openmoc/src/host/FlatSourceRegion.cpp',
                  'openmoc/src/host/Geometry.cpp',
                  'openmoc/src/host/LocalCoords.cpp',
                  'openmoc/src/host/log.cpp',
                  'openmoc/src/host/Material.cpp',
                  'openmoc/src/host/Point.cpp',
                  'openmoc/src/host/Quadrature.cpp',
                  'openmoc/src/host/Solver.cpp',
                  'openmoc/src/host/Surface.cpp',
                  'openmoc/src/host/Timer.cpp',
                  'openmoc/src/host/Track.cpp',
                  'openmoc/src/host/TrackGenerator.cpp',
                  'openmoc/src/host/Universe.cpp']

sources['cuda'] = ['openmoc/cuda/openmoc_cuda.i',
                   'openmoc/src/dev/DeviceMaterial.cu',
                   'openmoc/src/dev/DeviceTrack.cu',
                   'openmoc/src/dev/DeviceFlatSourceRegion.cu']



###############################################################################
#                                Compiler Paths
###############################################################################

path_to_gcc = '/usr/'
path_to_nvcc = '/usr/local/cuda-5.0/'
path_to_icpc = '/usr/intel/composer_xe_2013.1.117/composer_xe_2013.1.117/'

# Compiler binaries
gcc = path_to_gcc + 'bin/gcc'
icpc = path_to_icpc + 'bin/intel64/icpc'
nvcc = path_to_nvcc + 'bin/nvcc'



###############################################################################
#                                Compiler Flags
###############################################################################

compiler_flags = {}

compiler_flags['gcc'] = ['-c', 
                         '-O3', 
                         '-fopenmp', 
                         '-std=c++0x', 
                         '-fpic']

compiler_flags['icpc'] =['-c', 
                          '-O3', 
                          '-openmp', 
                          '-xhost', 
                          '-std=c++0x', 
                          '-fpic',
                          '-openmp-report', 
                          '-vec-report']

compiler_flags['nvcc'] =  ['-c', 
                           '-O3', 
                           '--ptxas-options=-v', 
                           '--compiler-options', 
                           '-fpic',
                           '-gencode=arch=compute_20,code=sm_20',
                           '-gencode=arch=compute_30,code=sm_30']



###############################################################################
#                                 Linker Flags
###############################################################################

linker_flags = {}

linker_flags['gcc'] = ['-lstdc++', 
                      '-lgomp', 
                      '-fopenmp', 
                      '-shared', 
                      '-Wl,-soname,_openmoc.so']

linker_flags['icpc'] = ['-lstdc++', 
                        '-openmp', 
                        '-liomp5', 
                        '-lpthread', 
                        '-lirc', 
                        '-limf', 
                        '-lrt', 
                        '-shared', 
                        '-Wl,-soname,_openmoc.so']

linker_flags['nvcc'] = ['-shared', 
                       'build/lib.linux-x86_64-2.7/_openmoc.so']



###############################################################################
#                               Shared Libraries
###############################################################################

shared_libraries = {}

shared_libraries['gcc'] = []
shared_libraries['icpc'] = []
shared_libraries['nvcc'] = ['cudart']



###############################################################################
#                              Library Directories
###############################################################################

library_directories = {}

library_directories['gcc'] = []
library_directories['icpc'] = [path_to_icpc + 'compiler/lib/intel64']
library_directories['nvcc'] = [path_to_nvcc + 'lib64']



###############################################################################
#                              Include Directories
###############################################################################

include_directories = {}

include_directories['gcc'] = []
include_directories['icpc'] =[path_to_icpc + 'compiler/include']
include_directories['nvcc'] = [path_to_nvcc + 'include']



###############################################################################
#                                 SWIG Flags
###############################################################################

swig_flags = ['-c++', '-keyword']



###############################################################################
#                                  Macros
###############################################################################

macros = {}
macros['gcc'] = {}
macros['icpc'] = {}
macros['nvcc'] = {}

macros['gcc']['single']= [('FP_PRECISION', 'float'), 
                          ('SINGLE', None),
                          ('GNU', None)]
macros['icpc']['single']= [('FP_PRECISION', 'float'), 
                           ('SINGLE', None),
                           ('INTEL', None)]

macros['gcc']['double'] = [('FP_PRECISION', 'double'), 
                           ('DOUBLE', None),
                           ('GNU', None)]
macros['icpc']['double'] = [('FP_PRECISION', 'double'), 
                            ('DOUBLE', None),
                            ('INTEL', None)]

macros['nvcc']['single'] = [('FP_PRECISION', 'float'), 
                            ('SINGLE', None),
                            ('CUDA', None)]
macros['nvcc']['double'] = [('FP_PRECISION', 'double'), 
                            ('DOUBLE', None),
                            ('CUDA', None)]
