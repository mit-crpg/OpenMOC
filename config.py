###############################################################################
#                                 User Options
###############################################################################

# Name of the Python package to be created
package_name = 'openmoc'

# Supported C++ compilers: 'gcc', 'icpc', 'all'
cpp_compilers = ['icpc']

default_cc = 'gcc'

# Supported floating point precision: 'single', 'double', 'all'
fp_precision = ['single']

default_fp = 'single'

# Compile using ccache (most relevant for developers)
with_ccache = True

# Compile with debug flags
debug_mode = True

# Use CUDA set to True or False
with_cuda = False

# Compile module for the Intel Xeon Phi (MIC)
with_mic = True



###############################################################################
#                                 Source Code
###############################################################################

sources = {}

sources['c++'] = ['openmoc/openmoc.i',
                  'src/Cell.cpp',
                  'src/Geometry.cpp',
                  'src/LocalCoords.cpp',
                  'src/log.cpp',
                  'src/Material.cpp',
                  'src/Point.cpp',
                  'src/Quadrature.cpp',
                  'src/Solver.cpp',
                  'src/CPUSolver.cpp',
                  'src/Surface.cpp',
                  'src/Timer.cpp',
                  'src/Track.cpp',
                  'src/TrackGenerator.cpp',
                  'src/Universe.cpp']

sources['cuda'] = ['openmoc/cuda/openmoc_cuda.i',
                   'src/dev/gpu/clone.cu',
                   'src/dev/gpu/GPUQuery.cu',
                   'src/dev/gpu/GPUSolver.cu']

sources['mic'] =  ['openmoc/mic/openmoc_mic.i',
                   'src/dev/mic/clone.cpp',
                   'src/dev/mic/MICQuery.cpp',
                   'src/dev/mic/MICSolver.cpp']




###############################################################################
#                                Compiler Paths
###############################################################################

path_to_gcc = '/usr/'
path_to_nvcc = '/usr/local/cuda/'
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
                           '-O3',
                           '--compiler-options', 
                           '-fpic',
                           '-gencode=arch=compute_20,code=sm_20',
                           '-gencode=arch=compute_30,code=sm_30',
                           '--ptxas-options=-v']



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
                        '-Xlinker',
                        '-soname=_openmoc.so']

linker_flags['mic'] = ['-lstdc++', 
                       '-openmp', 
                       '-liomp5', 
                       '-lpthread', 
                       '-lirc', 
                       '-limf', 
                       '-lrt',
                       '-shared',
                       'build/lib.linux-x86_64-2.7/_openmoc.so',
                       '-Xlinker',
                       '-soname=_openmoc.so',
                       '-loffload']

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
                            ('CUDA', None),
                            ('CCACHE_CC', 'nvcc')]
macros['nvcc']['double'] = [('FP_PRECISION', 'double'), 
                            ('DOUBLE', None),
                            ('CUDA', None),
                            ('CCACHE_CC', 'nvcc')]
