import sys, sysconfig
import copy
import numpy
from distutils.extension import Extension
from distutils.util import get_platform
from distutils.dist import Distribution
from distutils.command.install_lib import install_lib


def get_openmoc_object_name():
  """Returns the name of the main openmoc shared library object"""

  ext_suffix = sysconfig.get_config_var('SOABI')

  if ext_suffix is None:
    filename = '_openmoc.so'
  else:
    filename = '_openmoc.{0}.so'.format(ext_suffix)

  return filename


def get_shared_object_path():
  """Returns the name of the distutils build directory"""

  install_lib_command = install_lib(Distribution())
  install_lib_command.initialize_options()
  install_lib_command.finalize_options()

  directory = install_lib_command.build_dir

  return directory


def get_openmoc():
  """Returns the path and name of the main shared library object"""

  return get_shared_object_path() + '/' + get_openmoc_object_name()



class configuration:
  """User-defined build configuration options for OpenMOC

  Configuration options may be set using compile time flags. To view a
  list of these options, run 'python setup.py install --help' in the
  console. The default configuration options are shown below and should
  only be revised by developers familiar with the code and its configuration
  management system.
  """

  #############################################################################
  #                               User Options
  #############################################################################

  # Default C++ compiler for the main openmoc module is GCC
  cc = 'gcc'

  # Default floating point for the main openmoc module is single precision
  fp = 'single'

  # Compile using ccache (for developers needing fast recompilation)
  with_ccache = False

  # Compile code with debug symbols (ie, -g)
  debug_mode = False

  # Compile code with debug symbols (ie, -g, -pg)
  profile_mode = False

  # Build the openmoc.cuda module
  with_cuda = False

  # The vector length used for the VectorizedSolver class. This will used
  # as a hint for the Intel compiler to issue SIMD (ie, SSE, AVX, etc) vector
  # instructions. This is accomplished by adding "dummy" energy groups such
  # that the number of energy groups is be fit too a multiple of this
  # vector_length, and restructuring the innermost loops in the solver to
  # loop from 0 to the vector length
  vector_length = 8

  # The vector alignment used in the VectorizedSolver class when allocating
  # aligned data structures using MM_MALLOC and MM_FREE
  vector_alignment = 16

  # List of C/C++/CUDA distutils.extension objects which are created based
  # on which flags are specified at compile time.
  extensions = list()

  # List of the possible packages to install based on runtime options
  packages = ['openmoc', 'openmoc.cuda']


  #############################################################################
  #                                 Source Code
  #############################################################################

  # Dictionary of source code files to compile for each extension module
  sources = dict()

  sources['gcc'] = ['openmoc/openmoc_wrap.cpp',
                    'src/Cell.cpp',
                    'src/Geometry.cpp',
                    'src/LocalCoords.cpp',
                    'src/log.cpp',
                    'src/Material.cpp',
                    'src/MOCKernel.cpp',
                    'src/Point.cpp',
                    'src/Quadrature.cpp',
                    'src/ExpEvaluator.cpp',
                    'src/Solver.cpp',
                    'src/CPUSolver.cpp',
                    'src/Surface.cpp',
                    'src/Timer.cpp',
                    'src/Track.cpp',
                    'src/TrackGenerator.cpp',
                    'src/TrackTraversingAlgorithms.cpp',
                    'src/TraverseTracks.cpp',
                    'src/Universe.cpp',
                    'src/Vector.cpp',
                    'src/Matrix.cpp',
                    'src/Cmfd.cpp',
                    'src/linalg.cpp']

  sources['clang'] = ['openmoc/openmoc_wrap.cpp',
                      'src/Cell.cpp',
                      'src/Geometry.cpp',
                      'src/LocalCoords.cpp',
                      'src/log.cpp',
                      'src/Material.cpp',
                      'src/MOCKernel.cpp',
                      'src/Point.cpp',
                      'src/Quadrature.cpp',
                      'src/ExpEvaluator.cpp',
                      'src/Solver.cpp',
                      'src/CPUSolver.cpp',
                      'src/Surface.cpp',
                      'src/Timer.cpp',
                      'src/Track.cpp',
                      'src/TrackGenerator.cpp',
                      'src/TrackTraversingAlgorithms.cpp',
                      'src/TraverseTracks.cpp',
                      'src/Universe.cpp',
                      'src/Cmfd.cpp',
                      'src/Vector.cpp',
                      'src/Matrix.cpp',
                      'src/linalg.cpp']


  sources['icpc'] = ['openmoc/openmoc_wrap.cpp',
                     'src/Cell.cpp',
                     'src/Geometry.cpp',
                     'src/LocalCoords.cpp',
                     'src/log.cpp',
                     'src/Material.cpp',
                     'src/MOCKernel.cpp',
                     'src/Point.cpp',
                     'src/Quadrature.cpp',
                     'src/ExpEvaluator.cpp',
                     'src/Solver.cpp',
                     'src/CPUSolver.cpp',
                     'src/VectorizedSolver.cpp',
                     'src/Surface.cpp',
                     'src/Timer.cpp',
                     'src/Track.cpp',
                     'src/TrackGenerator.cpp',
                     'src/TrackTraversingAlgorithms.cpp',
                     'src/TraverseTracks.cpp',
                     'src/Universe.cpp',
                     'src/Cmfd.cpp',
                     'src/Vector.cpp',
                     'src/Matrix.cpp',
                     'src/linalg.cpp']


  sources['bgxlc'] = ['openmoc/openmoc_wrap.cpp',
                      'src/Cell.cpp',
                      'src/Geometry.cpp',
                      'src/LocalCoords.cpp',
                      'src/log.cpp',
                      'src/Material.cpp',
                      'src/MOCKernel.cpp',
                      'src/Point.cpp',
                      'src/Quadrature.cpp',
                      'src/ExpEvaluator.cpp',
                      'src/Solver.cpp',
                      'src/CPUSolver.cpp',
                      'src/Surface.cpp',
                      'src/Timer.cpp',
                      'src/Track.cpp',
                      'src/TrackGenerator.cpp',
                      'src/TrackTraversingAlgorithms.cpp',
                      'src/TraverseTracks.cpp',
                      'src/Universe.cpp',
                      'src/Cmfd.cpp',
                      'src/Vector.cpp',
                      'src/Matrix.cpp',
                      'src/linalg.cpp']


  sources['nvcc'] = ['openmoc/cuda/openmoc_cuda_wrap.cpp',
                     'src/accel/cuda/GPUExpEvaluator.cu',
                     'src/accel/cuda/GPUQuery.cu',
                     'src/accel/cuda/clone.cu',
                     'src/accel/cuda/GPUSolver.cu']


  #############################################################################
  #                                Compiler Flags
  #############################################################################

  # A dictionary of the compiler flags to use for each compiler type
  compiler_flags = dict()

  compiler_flags['gcc'] = ['-c', '-O3', '-ffast-math', '-fopenmp',
                           '-std=c++11', '-fpic']
  compiler_flags['clang'] = ['-c', '-O3', '-ffast-math', '-std=c++11',
                             '-fopenmp', '-fvectorize', '-fpic',
                             '-Qunused-arguments',
                             '-Wno-deprecated-register',
                             '-Wno-parentheses-equality']
  compiler_flags['icpc'] =['-c', '-O3', '-fast', '--ccache-skip',
                           '-openmp', '-xhost', '-std=c++11',
                           '--ccache-skip', '-fpic',
                           '-openmp-report', '-vec-report']
  compiler_flags['bgxlc'] = ['-c', '-O2', '-qarch=qp', '-qreport',
                             '-qsimd=auto', '-qtune=qp', '-qunroll=auto',
                             '-qsmp=omp', '-qpic']
  compiler_flags['nvcc'] =  ['--relocatable-device-code', 'true',
                             '-c', '-O3',  '-std=c++11',
                             '--compiler-options', '-fpic',
                             '-arch=compute_20']


  #############################################################################
  #                                 Linker Flags
  #############################################################################

  # A dictionary of the linker flags to use for each compiler type
  linker_flags = dict()

  if ('macosx' in get_platform()):
    linker_flags['gcc'] = ['-fopenmp', '-dynamiclib', '-lpython2.7',
                           '-Wl,-install_name,' + get_openmoc_object_name()]
  else:
    linker_flags['gcc'] = ['-fopenmp', '-shared',
                           '-Wl,-soname,' + get_openmoc_object_name()]

  if ('macosx' in get_platform()):
    linker_flags['clang'] = ['-fopenmp', '-dynamiclib', '-lpython2.7',
                             '-Wl,-install_name,' + get_openmoc_object_name()]
  else:
    linker_flags['clang'] = ['-fopenmp', '-shared',
                             '-Wl,-soname,' + get_openmoc_object_name()]

  linker_flags['icpc'] = [ '-openmp', '-shared',
                           '-Xlinker', '-soname=' + get_openmoc_object_name()]
  linker_flags['bgxlc'] = ['-qmkshrobj', '-shared',
                           '-R/soft/compilers/ibmcmp-may2013/lib64/bg/bglib64',
                           '-Wl,-soname,' + get_openmoc_object_name()]
  linker_flags['nvcc'] = ['-shared', get_openmoc()]


  #############################################################################
  #                               Shared Libraries
  #############################################################################

  # A dictionary of the shared libraries to use for each compiler type
  shared_libraries = dict()

  shared_libraries['gcc'] = ['stdc++', 'gomp', 'dl','pthread', 'm']
  shared_libraries['clang'] = ['stdc++', 'gomp', 'dl','pthread', 'm']
  shared_libraries['icpc'] = ['stdc++', 'iomp5', 'pthread', 'irc',
                              'imf','rt', 'mkl_rt','m',]
  shared_libraries['bgxlc'] = ['stdc++', 'pthread', 'm', 'xlsmp', 'rt']
  shared_libraries['nvcc'] = ['cudadevrt', 'cudart']


  #############################################################################
  #                              Library Directories
  #############################################################################

  # A dictionary of the library directories to use for each compiler type
  # if not set in the LD_LIBRARY_PATH environment variable
  library_directories = dict()

  usr_lib = sys.exec_prefix + '/lib'

  library_directories['gcc'] = [usr_lib]
  library_directories['clang'] = [usr_lib]
  library_directories['icpc'] = [usr_lib]
  library_directories['bgxlc'] = [usr_lib]
  library_directories['nvcc'] = [usr_lib, '/usr/local/cuda/lib64']


  #############################################################################
  #                              Include Directories
  #############################################################################

  # A dictionary of the include directories to use for each compiler type
  # for header files not found from paths set in the user's environment
  include_directories = dict()

  include_directories['gcc'] = list()
  include_directories['clang'] = list()
  include_directories['icpc'] = list()
  include_directories['bgxlc'] = list()
  include_directories['nvcc'] = ['/usr/local/cuda/include']


  ###########################################################################
  #                                 SWIG Flags
  ###########################################################################

  # A list of the flags for SWIG
  swig_flags = ['-c++', '-python', '-keyword']

  # Python 3 only
  if sys.version_info[0] == 3:
    swig_flags.append('-py3')


  #############################################################################
  #                                  Macros
  #############################################################################

  # A dictionary of the macros to set at compile time for each compiler type
  # and floating point precisin level
  macros = dict()

  macros['gcc'] = dict()
  macros['clang'] = dict()
  macros['icpc'] = dict()
  macros['bgxlc'] = dict()
  macros['nvcc'] = dict()

  macros['gcc']['single']= [('FP_PRECISION', 'float'),
                            ('SINGLE', None),
                            ('GCC', None),
                            ('VEC_LENGTH', vector_length),
                            ('VEC_ALIGNMENT', vector_alignment)]

  macros['clang']['single']= [('FP_PRECISION', 'float'),
                              ('SINGLE', None),
                              ('CLANG', None),
                              ('VEC_LENGTH', vector_length),
                              ('VEC_ALIGNMENT', vector_alignment)]

  macros['icpc']['single']= [('FP_PRECISION', 'float'),
                             ('SINGLE', None),
                             ('ICPC', None),
                             ('MKL_ILP64', None),
                             ('VEC_LENGTH', vector_length),
                             ('VEC_ALIGNMENT', vector_alignment)]

  macros['bgxlc']['single'] = [('FP_PRECISION', 'float'),
                               ('SINGLE', None),
                               ('BGXLC', None),
                               ('VEC_LENGTH', vector_length),
                               ('VEC_ALIGNMENT', vector_alignment),
                               ('CCACHE_CC', 'bgxlc++_r')]

  macros['nvcc']['single'] = [('FP_PRECISION', 'float'),
                              ('SINGLE', None),
                              ('NVCC', None),
                              ('CCACHE_CC', 'nvcc')]

  macros['gcc']['double'] = [('FP_PRECISION', 'double'),
                             ('DOUBLE', None),
                             ('GCC', None),
                             ('VEC_LENGTH', vector_length),
                             ('VEC_ALIGNMENT', vector_alignment)]

  macros['clang']['double'] = [('FP_PRECISION', 'double'),
                               ('DOUBLE', None),
                               ('CLANG', None),
                               ('VEC_LENGTH', vector_length),
                               ('VEC_ALIGNMENT', vector_alignment)]

  macros['icpc']['double'] = [('FP_PRECISION', 'double'),
                              ('DOUBLE', None),
                              ('ICPC', None),
                              ('MKL_ILP64', None),
                              ('VEC_LENGTH', vector_length),
                              ('VEC_ALIGNMENT', vector_alignment)]

  macros['bgxlc']['double'] = [('FP_PRECISION', 'double'),
                               ('DOUBLE', None),
                               ('BGXLC', None),
                               ('VEC_LENGTH', vector_length),
                               ('VEC_ALIGNMENT', vector_alignment),
                               ('CCACHE_CC', 'bgxlc++_r')]

  macros['nvcc']['double'] = [('FP_PRECISION', 'double'),
                              ('DOUBLE', None),
                              ('NVCC', None),
                              ('CCACHE_CC', 'nvcc')]

  # define OPENMP and SWIG (for log output)
  for compiler in macros:
    for precision in macros[compiler]:
      macros[compiler][precision].append(('OPENMP', None))
      macros[compiler][precision].append(('SWIG', None))


  def setup_extension_modules(self):
    """Sets up the C/C++/CUDA extension modules for this distribution.

    Create list of extensions for Python modules within the openmoc
    Python package based on the user-defined flags defined at compile time.
    """

    # If the user wishes to compile using debug mode, append the debugging
    # flag to all lists of compiler flags for all distribution types
    if self.debug_mode:
      for k in self.compiler_flags:
        self.compiler_flags[k].append('-g')

    # If the user wishes to compile using profile mode, append the profiling
    # flag to all lists of compiler flags for all distribution types
    if self.profile_mode:
      for k in self.compiler_flags:
        self.compiler_flags[k].append('-pg')
        self.compiler_flags[k].append('-g')

    # Obtain the NumPy include directory
    try:
      numpy_include = numpy.get_include()
    except AttributeError:
      numpy_include = numpy.get_numpy_include()

    # Add the NumPy include directory to the include directories
    # list for each type of compiler
    for cc in self.include_directories.keys():
      self.include_directories[cc].append(numpy_include)


    # The main openmoc extension (defaults are gcc and single precision)
    self.swig_flags += ['-D' + self.fp.upper()]
    if self.fp == 'double':
      self.swig_flags += ['-DFP_PRECISION=double']
    else:
      self.swig_flags += ['-DFP_PRECISION=float']

    self.extensions.append(
      Extension(name = '_openmoc',
                sources = copy.deepcopy(self.sources[self.cc]),
                library_dirs = self.library_directories[self.cc],
                libraries = self.shared_libraries[self.cc],
                extra_link_args = self.linker_flags[self.cc],
                include_dirs = self.include_directories[self.cc],
                define_macros = self.macros[self.cc][self.fp],
                swig_opts = self.swig_flags + ['-D' + self.cc.upper()]))

    # The openmoc.cuda extension if requested by the user at compile
    # time (--with-cuda)
    if self.with_cuda:

      self.extensions.append(
        Extension(name = '_openmoc_cuda',
                  sources = copy.deepcopy(self.sources['nvcc']),
                  library_dirs = self.library_directories['nvcc'],
                  libraries = self.shared_libraries['nvcc'],
                  extra_link_args = self.linker_flags['nvcc'],
                  include_dirs = self.include_directories['nvcc'],
                  define_macros = self.macros['nvcc'][self.fp],
                  swig_opts = self.swig_flags  + ['-DNVCC'],
                  export_symbols = ['init_openmoc']))
