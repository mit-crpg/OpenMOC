##############################################################################
#                                 User Options
###############################################################################

# Name of the Python package to be created
package_name = 'openmoc'

# Supported C++ compilers: 'gcc', 'icpc', 'all'
cpp_compilers = ['icpc']

# Supported floating point precision: 'single', 'double', 'all'
fp_precision = ['single']

# Default floating point precision to use for main openmoc module
default_fp = 'double'

# Compile using ccache (most relevant for developers)
with_ccache = True

# Compile with debug flags
debug_mode = False

# Use CUDA set to True or False
with_cuda = False

# Compile package for Intel MIC
with_mic = True



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
                   'openmoc/src/dev/DeviceTrack.cu',
                   'openmoc/src/dev/DeviceMaterial.cu',
                   'openmoc/src/dev/DeviceQuery.cu',
                   'openmoc/src/dev/DeviceSolver.cu']

sources['mic'] = ['openmoc/mic/openmoc_mic.i',
                  'openmoc/src/mic/DeviceQuery.cpp']



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
                          '--ccache-skip',
                          '-openmp', 
                          '--ccache-skip',
                          '-xhost', 
                          '-std=c++0x', 
                          '-fpic',
                          '--ccache-skip',
                          '-openmp-report', 
                          '-vec-report']

compiler_flags['nvcc'] =  ['-c', 
                           '-O0',
                           '--ccache-skip',
                           '--compiler-options', 
                           '-fpic',
                           '--ccache-skip',
                           '-gencode=arch=compute_20,code=sm_20',
                           '--ccache-skip',
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
                        'build/lib.linux-x86_64-2.6/_openmoc.so',
                        '-Xlinker',
                        '-soname=_openmoc.so',
                        '-loffload']

linker_flags['nvcc'] = ['-shared', 
                       'build/lib.linux-x86_64-2.6/_openmoc.so']



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
library_directories['icpc'] = []
library_directories['nvcc'] = []




###############################################################################
#                              Include Directories
###############################################################################

include_directories = {}

include_directories['gcc'] = []
include_directories['icpc'] = []
include_directories['nvcc'] = []



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
                            ('CUDA', None),
                            ('CCACHE_CC', 'nvcc')]
macros['nvcc']['double'] = [('FP_PRECISION', 'double'), 
                            ('DOUBLE', None),
                            ('CUDA', None),
                            ('CCACHE_CC', 'nvcc')]
