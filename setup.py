'''
The setup script for OpenMOC
'''

import os, shutil, site, string, sys
from distutils.errors import DistutilsOptionError
import distutils.ccompiler
import multiprocessing
import multiprocessing.pool as pool

# Use setuptools only if the user opts-in by setting the USE_SETUPTOOLS.
# This ensures consistent behavior but allows for advanced usage with
# virtualenv, buildout, pip and others.
with_setuptools = False
if 'USE_SETUPTOOLS' in os.environ or 'pip' in __file__:
  with_setuptools = True

if with_setuptools:
  from setuptools import setup
  from setuptools.command.build_ext import build_ext
  from setuptools.command.build_py import build_py
  from setuptools.command.install import install
else:
  from distutils.core import setup
  from distutils.command.build_ext import build_ext
  from distutils.command.build_py import build_py
  from distutils.command.install import install

import config


# Instantiate a configurations class for this OpenMOC build
config = config.configuration()


class custom_install(install):
  """Defines the compile time options for OpenMOC.

  This class derives from the distutils.command.install class. Distutils
  provides a set of flags which may be invoked by the user at compile
  time. The custom_install class adds to that list a series of options
  which are useful in compiling a specific distribution of OpenMOC.

  To view a list of available compile time options, simply type the following
  from a console:

  python setup.py install --help

  The custom_install class extracts user-defined options from the command
  line and uses them to populate the config.configurations class. The
  configurations class then uses these options to generate the Python
  C/C++/CUDA extension objects which are delivered to the distutils
  setup method at the end of this script.

  Developers wishing to extend OpenMOC's functionality with new compilation
  options will need to extend this class and the configurations class.
  """

  # The user options for a customized OpenMOC build
  user_options = [
    ('cc=', None, "Compiler (gcc, icpc, bgxlc, mpicc) for main openmoc module"),
    ('fp=', None, "Floating point precision (single or double) for " + \
                  "main openmoc module"),
    ('ng=', None, "Specify number of groups (optional only for optimization)" + \
                  " for main openmoc module"),
    ('with-cuda', None, "Build openmoc.cuda module for NVIDIA GPUs"),
    ('debug-mode', None, "Build with debugging symbols"),
    ('sanitizer-mode', None, "Build with address sanitizer"),
    ('profile-mode', None, "Build with profiling symbols"),
    ('coverage-mode', None, "Build with coverage symbols"),
    ('with-ccache', None, "Build with ccache for rapid recompilation"),
  ]

  # Include all of the default options provided by distutils for the
  # install command parent class
  user_options += install.user_options

  # Set some compile options to be boolean switches
  boolean_options = ['debug-mode',
                     'sanitizer-mode',
                     'profile-mode',
                     'coverage-mode',
                     'with-ccache']

  # Include all of the boolean options provided by distutils for the
  # install command parent class
  boolean_options += install.boolean_options


  def initialize_options(self):
    """Set the default OpenMOC build options

    The default installation is invoked by following console command:

    python setup.py install

    By default, this will build the main openmoc C/C++ Python extension using
    the GCC compiler with single precision (default).
    """

    # Run the install command parent class' initialize_options method
    install.initialize_options(self)

    # Default compiler and precision level for the main openmoc module
    self.cc = 'gcc'
    self.fp = 'single'
    self.ng = None

    # Set defaults for each of the newly defined compile time options
    self.with_cuda = False
    self.debug_mode = False
    self.sanitizer_mode = False
    self.profile_mode = False
    self.coverage_mode = False
    self.with_ccache = False

    # Check that swig is installed:
    swigLocation = shutil.which('swig')
    if not swigLocation:
        print("-> Please install swig before building <-")
        sys.exit()


  def finalize_options(self):
    """Extract options from the flags invoked by the user at compile time.

    This method performs error checking of the options specified by
    the user at compile time, and initialize the config.configurations
    class instance. The method conclude with a call to the
    configurations.setup_extension_modules class method which creates
    the C/C++/CUDA extension modules to be passed to the distutils
    setup method at the end of this script.
    """

    # Run the install command parent class' finalize_options method
    install.finalize_options(self)

    # Set the configuration options specified to be the default
    # unless the corresponding flag was invoked by the user
    config.num_groups = self.ng
    config.with_cuda = self.with_cuda
    config.debug_mode = self.debug_mode
    config.sanitizer_mode = self.sanitizer_mode
    config.profile_mode = self.profile_mode
    config.coverage_mode = self.coverage_mode
    config.with_ccache = self.with_ccache

    # Check that the user specified a supported C++ compiler
    if self.cc not in ['gcc', 'clang', 'icpc', 'bgxlc', 'mpicc']:
      raise DistutilsOptionError \
            ('Must supply the -cc flag with one of the supported ' +
             'C++ compilers: gcc, clang, icpc, bgxlc, mpicc')
    else:
      config.cc = self.cc

    # Check that the user specified a supported floating point precision
    if self.fp not in ['single', 'double']:
      raise DistutilsOptionError \
          ('Must supply the -cc flag with one of the supported ' +
           'floating point precision levels: single, double')
    else:
      config.fp = self.fp

    # Build the C/C++/CUDA extension modules for this distribution
    config.setup_extension_modules()



def customize_compiler(self):
  """Inject redefined _compile method into distutils

  This method enables us to choose compilers based on the macros defined
  in the compiler flags (ie, '-DGCC', '-DNVCC', etc), or on the
  source extension (ie, *.cpp, *.cu, etc.).

  Adapted from Robert McGibbon's CUDA distutils setup provided in open source
  form here: https://github.com/rmcgibbo/npcuda-example
  """

  # Inform the compiler it can processes .cu CUDA source files
  self.src_extensions.append('.cu')

  # Save reference to the default _compile method
  super_compile = self._compile

  # Redefine the _compile method. This gets executed for each
  # object but distutils doesn't have the ability to change compilers
  # based on source extension, so we add that functionality here
  def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):

    # If compiler is GNU's gcc and the source is C++, use gcc
    if config.cc == 'gcc' and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache gcc')
      else:
        self.set_executable('compiler_so', 'gcc')

      postargs = config.compiler_flags['gcc']

    # If compiler is GNU's gcc and the source is C++, use gcc
    elif config.cc == 'mpicc' and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache mpicc')
      else:
        self.set_executable('compiler_so', 'mpicc')

      postargs = config.compiler_flags['mpicc']



    # If compiler is Apple's clang and the source is C++, use clang
    elif config.cc == 'clang' and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache clang')
      else:
        self.set_executable('compiler_so', 'clang')

      postargs = config.compiler_flags['clang']

    # If compiler is Intel's icpc and the source is C++, use icpc
    elif config.cc == 'icpc' and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache icpc')
      else:
        self.set_executable('compiler_so', 'icpc')

      postargs = config.compiler_flags['icpc']

    # If compiler is IBM's bgxlc and the source is C++, use bgxlc
    elif config.cc == 'bgxlc' and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'bgxlc++_r')
      else:
        self.set_executable('compiler_so', 'bgxlc++_r')

      postargs = config.compiler_flags['bgxlc']

    # If NVCC is a defined macro and the source is C++, compile
    # SWIG-wrapped CUDA code with gcc
    elif '-DNVCC' in pp_opts and os.path.splitext(src)[1] == '.cpp':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache gcc')
      else:
        self.set_executable('compiler_so', 'gcc')

      postargs = config.compiler_flags['gcc']

    # If NVCC is a defined macro and the source is CUDA, use nvcc
    elif '-DNVCC' in pp_opts and os.path.splitext(src)[1] == '.cu':
      if config.with_ccache:
        self.set_executable('compiler_so', 'ccache nvcc')
      else:
        self.set_executable('compiler_so', 'nvcc')

      postargs = config.compiler_flags['nvcc']

    # If we cannot determine how to compile this file, throw exception
    else:
      raise EnvironmentError('Unable to compile ' + str(src))

    # Now call distutils-defined _compile method
    super_compile(obj, src, ext, cc_args, postargs, pp_opts)

  # Inject our redefined _compile method into the class
  self._compile = _compile


def customize_linker(self):
  """Inject redefined link method into distutils

  This method enables us to choose the linker based on the name
  of the output shared library filename (ie, _openmoc_intel_single.so)

  Adapted from Robert McGibbon's CUDA distutils setup provided in open source
  form here: https://github.com/rmcgibbo/npcuda-example
  """

  # Save references to the default link method
  super_link = self.link

  # Redefine the link method. This gets executed to link each extension
  # module. We add the functionality to choose the compiler for linking
  # based on the name of the extension module
  def link(target_desc, objects, output_filename,
           output_dir=None, libraries=None,
           library_dirs=None, runtime_library_dirs=None,
           export_symbols=None, debug=0, extra_preargs=None,
           extra_postargs=None, build_temp=None, target_lang=None):

    if config.cc == 'gcc':
      self.set_executable('linker_so', 'gcc')
      self.set_executable('linker_exe', 'gcc')

    elif config.cc == 'mpicc':
      self.set_executable('linker_so', 'mpicc')
      self.set_executable('linker_exe', 'mpicc')

    elif config.cc == 'clang':
      self.set_executable('linker_so', 'clang')
      self.set_executable('linker_exe', 'clang')

    elif config.cc == 'icpc':
      self.set_executable('linker_so', 'icpc')
      self.set_executable('linker_exe', 'icpc')

    elif config.cc == 'bgxlc':
      self.set_executable('linker_so', 'bgxlc++_r')
      self.set_executable('linker_exe', 'bgxlc++_r')

    # If the filename for the extension contains cuda, use nvcc to link
    if 'cuda' in output_filename:
      self.set_executable('linker_so', 'nvcc')
      self.set_executable('linker_exe', 'nvcc')

    # Now call distutils-defined link method
    super_link(target_desc, objects,
               output_filename, output_dir, libraries,
               library_dirs, runtime_library_dirs,
               export_symbols, debug, extra_preargs,
               extra_postargs, build_temp)

  # Inject our redefined link method into the class
  self.link = link


# monkey-patch for parallel compilation
def parallel_compile(self, sources, output_dir=None, macros=None,
                     include_dirs=None, debug=0, extra_preargs=None,
                     extra_postargs=None, depends=None):
  """A parallel version of the Distutils compile method

  Note that this routine is modified from StackOverflow post #11013851
  """

  # Copy args from distutils.ccompiler.CCompiler directly
  macros, objects, extra_postargs, pp_opts, build = \
       self._setup_compile(output_dir, macros, include_dirs,
                           sources, depends, extra_postargs)
  cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
  num_cpus = multiprocessing.cpu_count()

  # Define routine for each thread to use to compile on its own
  def _single_compile(obj):
      try: src, ext = build[obj]
      except KeyError: return
      self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

  # Convert thread mapping to C/C++/CUDA objects to a list and return
  list(pool.ThreadPool(num_cpus).map(_single_compile, objects))
  return objects


# Inject parallel_compile to override distutils sequential compile method
distutils.ccompiler.CCompiler.compile=parallel_compile


# Run the customize_compiler to inject redefined and customized _compile and
# link methods into distutils
class custom_build_ext(build_ext):
  """Customizes distutils to work with different compiler types

  This class derives from the distutils.command.build_ext command class.
  It extends build_ex by creates customized compile and link methods
  which can accommodate different compiler types and options.
  """

  def build_extensions(self):
    customize_compiler(self.compiler)
    customize_linker(self.compiler)

    # Append a macro for the compiler to the SWIG flags
    swig_flags = config.swig_flags + ['-D' + config.cc.upper()]

    os.system('swig {0} -o '.format(str.join(' ', swig_flags)) + \
              'openmoc/swig/openmoc_wrap.cpp openmoc/swig/openmoc.i')

    if config.with_cuda:
      swig_flags = config.swig_flags + ['-DNVCC']
      os.system('swig {0} -o '.format(str.join(' ', swig_flags)) + \
                'openmoc/cuda/openmoc_cuda_wrap.cpp ' + \
                'openmoc/cuda/openmoc_cuda.i')

    # Move openmoc.py file created by swig into main python API folder
    os.system('mv openmoc/swig/openmoc.py openmoc/openmoc.py')

    build_ext.build_extensions(self)


# Run the distutils setup method for the complete build
dist = setup(name = 'openmoc',
      version = '0.4.0',
      description = 'An open source method of characteristics code for ' + \
                    'solving the 2D neutron distribution in nuclear reactors',
      url = 'https://github.com/mit-crpg/OpenMOC',
      download_url = 'https://github.com/mit-crpg/OpenMOC/tarball/v0.4.0',

      # Set the C/C++/CUDA extension modules built in setup_extension_modules()
      # in config.py based on the user-defined flags at compile time
      ext_modules = config.extensions,

      # Extract all of the Python packages for OpenMOC
      # (ie, openmoc.log, openmoc.materialize, etc)
      packages = config.packages,

      # Include SWIG interface files in package - this is important for PyPI
      package_data = {'' : ['*.i*']},

      # Inject our custom compiler and linker triggers
      cmdclass={ 'build_ext': custom_build_ext,
                 'install': custom_install}
)

# Rerun the build_py to setup links for C++ extension modules created by SWIG
# This prevents us from having to install twice
build_py = build_py(dist)
build_py.ensure_finalized()
build_py.run()

# Remove the shared library in the site packages
if "clean" in sys.argv:
    install_location = site.getsitepackages()[0]
    print("Removing build from "+ install_location)
    os.system("rm -rf " + install_location + "/*openmoc*")
    install_location = site.getusersitepackages()
    print("Removing build from "+ install_location)
    os.system("rm -rf " + install_location + "/*openmoc*")
    install_location = "./tests/"
    print("Removing build from "+ install_location)
    os.system("rm -rf "+install_location+"openmoc "+install_location+"build")
