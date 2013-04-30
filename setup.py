import os
from os.path import join as pjoin
from setuptools import setup, find_packages
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
import subprocess
import numpy

# Adapted from Robert McGibbon's CUDA distutils setup provide in open source 
# form here: https://github.com/rmcgibbo/npcuda-example


cpp_compiler = 'all'
fp_precision = 'all'
use_cuda = False


def find_in_path(name, path):
    "Find a file in a search path"
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_gcc():

    # Search the PATH for g++
    gcc = find_in_path('gcc', os.environ['PATH'])

    if gcc is None:
        raise EnvironmentError('The gcc binary could not be located in ' + \
                               'your $PATH. Please add it to your path.')

    gcc_home = os.path.dirname(os.path.dirname(os.path.dirname(gcc)))

    gccconfig = {'gcc':gcc}

    for k, v in gccconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The GNU %s path could not be located ' + \
                                   'in %s' % (k, v))

    return gccconfig


def locate_icpc():

    # Search the PATH for icpc
    icpc = find_in_path('icpc', os.environ['PATH'])

    if icpc is None:
        raise EnvironmentError('The icpc binary could not be located in ' + \
                               'your $PATH. Please add it to your path.')

    icpc_home = os.path.dirname(os.path.dirname(os.path.dirname(icpc)))

    icpcconfig = {'icpc':icpc, 'include': pjoin(icpc_home, 'compiler/include'),
                  'lib64': pjoin(icpc_home, 'compiler/lib/intel64')}

    for k, v in icpcconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The Intel %s path could not be located ' + \
                                   'in %s' % (k, v))

    return icpcconfig


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # Search the PATH for nvcc
    nvcc = find_in_path('nvcc', os.environ['PATH'])

    if nvcc is None:
        raise EnvironmentError('The nvcc binary could not be located in ' + \
                               'your $PATH. Please add it to your path.')
        
    cuda_home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'nvcc':nvcc, 'include': pjoin(cuda_home, 'include'),
                  'lib64': pjoin(cuda_home, 'lib64')}

    for k, v in cudaconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located ' + \
                                   'in %s' % (k, v))

    return cudaconfig


GCC  = locate_gcc()
ICPC = locate_icpc()
CUDA = locate_cuda()


# Obtain the numpy include directory
try:
    numpy_include = numpy.get_include()

except AttributeError:
    numpy_include = numpy.get_numpy_include()


sources = {}
sources['regular'] = ['openmoc/openmoc.i',
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

sources['cuda'] = ['openmoc/openmoc.i',
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
                   'openmoc/src/host/Universe.cpp',
                   'openmoc/src/dev/DeviceMaterial.cu',
                   'openmoc/src/dev/DeviceTrack.cu',
                   'openmoc/src/dev/DeviceFlatSourceRegion.cu']

swig_interface_file = 'openmoc'
swig_pkg_name = '_openmoc'
pkg_name = 'openmoc'
distro = ''
swig_opts = ['-c++', '-keyword']

extra_link_args = {}
extra_compile_args = {}
libraries = {}
library_dirs = {}
include_dirs = {}
runtime_library_dirs = {}
macros = {}

extra_link_args['gnu'] = ['-lstdc++', '-lgomp', '-fopenmp', '-shared',
                          '-Wl,-Bsymbolic-functions', '-Wl,-z,relro']
extra_link_args['intel'] = ['-lstdc++', '-openmp', '-liomp5', '-lpthread', 
                            '-lirc', '-limf', '-lrt', '-shared',
                            '-Wl,-Bsymbolic-functions', '-Wl,-z,relro']
extra_compile_args['gnu'] = ['-c', '-O3', '-fopenmp', '-std=c++0x', '-fPIC']
extra_compile_args['icpc'] =['-c', '-O3', '-openmp', '-std=c++0x', '-fpic',
                             '-xhost', '-openmp-report', '-vec-report']

libraries['gnu'] = []
libraries['intel'] = []

library_dirs['gnu'] = []
library_dirs['intel'] = [ICPC['lib64']]

include_dirs['gnu'] = [numpy_include, 'openmoc/src/host', 'openmoc/src/dev']
include_dirs['intel'] =[numpy_include, 'openmoc/src/host', 'openmoc/src/dev', 
                        ICPC['lib64']]

runtime_library_dirs['gnu'] = []
runtime_library_dirs['intel'] = [ICPC['lib64']]

macros = {}
macros['gnu'] = {}
macros['intel'] = {}
macros['gnu']['single']= [('FP_PRECISION', 'float'), 
                          ('SINGLE', None),
                          ('GNU', None)]
macros['intel']['single']= [('FP_PRECISION', 'float'), 
                            ('SINGLE', None),
                            ('INTEL', None)]
macros['gnu']['double'] = [('FP_PRECISION', 'double'), 
                           ('DOUBLE', None),
                           ('GNU', None)]
macros['intel']['double'] = [('FP_PRECISION', 'double'), 
                           ('DOUBLE', None),
                           ('INTEL', None)]


if cpp_compiler == 'gnu':
    cpp_compiler = [cpp_compiler]
    distro += '_gnu'

elif cpp_compiler == 'intel':
    cpp_compiler = [cpp_compiler]
    distro += '_intel'

elif cpp_compiler == 'all':
    cpp_compiler = ['gnu', 'intel']


if use_cuda:
    libraries['gnu'].extend(['cudart'])
    libraries['intel'].extend(['cudart'])
    library_dirs['gnu'].extend([CUDA['lib64']])
    library_dirs['intel'].extend([CUDA['lib64']])
    runtime_library_dirs['gnu'].extend([CUDA['lib64']])
    runtime_library_dirs['intel'].extend([CUDA['lib64']])
    extra_compile_args['nvcc'] = ['-c', '-arch=sm_20', '--ptxas-options=-v', 
                                  '--compiler-options', '-fpic']
    include_dirs['gnu'].extend([CUDA['include']])
    include_dirs['intel'].extend([CUDA['include']])
    sources = sources['cuda']

    distro += '_cuda'

else:
    sources = sources['regular']

# Define macros for single or double precision OpenMOC solver
if fp_precision == 'single':
    fp_precision = [fp_precision]
    distro += '_single'

elif fp_precision == 'double':
    fp_precision = [fp_precision]
    distro += '_double'

elif fp_precision == 'all':
    fp_precision = ['double', 'single']

swig_pkg_name += distro
pkg_name += distro
swig_interface_file += distro + '.i'

#sources.extend(['openmoc/' + swig_interface_file])

extensions = []
extensions.append(Extension(name = '_openmoc', 
                    sources = sources, 
                    library_dirs = library_dirs['gnu'], 
                    libraries = libraries['gnu'],
                    runtime_library_dirs = runtime_library_dirs['gnu'],
                    extra_link_args = extra_link_args['gnu'], 
                    include_dirs = include_dirs['gnu'],
                    define_macros = macros['gnu']['double'],
                    swig_opts = swig_opts))

sources.remove('openmoc/openmoc.i')
distro = ''

for fp in fp_precision:
    for cc in cpp_compiler:
        print fp
        print cc
        distro = 'openmoc_' + cc + '_' + fp 
        swig_interface_file = distro.replace('_', '/') + '/' + distro + '.i'
        sources.append(swig_interface_file)
        ext_name = '_' + distro.replace('.', '_')
#        ext_name = distro

        print distro
        print swig_interface_file

        extensions.append(Extension(name = ext_name, 
                    sources = sources, 
                    library_dirs = library_dirs[cc], 
                    libraries = libraries[cc],
                    runtime_library_dirs = runtime_library_dirs[cc],
                    extra_link_args = extra_link_args[cc], 
                    include_dirs = include_dirs[cc],
                    define_macros = macros[cc][fp],
                    swig_opts = swig_opts))

        distro = ''
#        sources.remove(swig_interface_file)


#openmoc = Extension(name = swig_pkg_name, 
#                    sources = sources, 
#                    library_dirs = library_dirs, 
#                    libraries = libraries,
#                    runtime_library_dirs = runtime_library_dirs,
#                    extra_link_args = extra_link_args, 
#                    include_dirs = include_dirs,
#                    define_macros = define_macros,
#                    swig_opts = swig_opts)


def customize_compiler(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""
    
    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super_compile = self._compile
    super_link = self.link

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, openmoc, cc_args, extra_postargs, pp_opts):

        print pp_opts
        print obj

        if '-DGNU' in pp_opts and os.path.splitext(src)[1] == '.cpp':
            self.set_executable('compiler_so', 'gcc')
            postargs = extra_compile_args['gnu']

        elif '-DINTEL' in pp_opts and os.path.splitext(src)[1] == '.cpp':
            self.set_executable('compiler_so', 'icpc')
            postargs = extra_compile_args['icpc']

        elif os.path.splitext(src)[1] == '.cu':
            self.set_executable('compiler_so', 'nvcc')
            postargs = extra_compile_args['nvcc']

            for item in extra_postargs:
                if 'intel' in item:
                    cc_args.remove(item)

        else:
            raise EnvironmentError('Unable to compile ' + str(src))

        super_compile(obj, src, openmoc, cc_args, postargs, pp_opts)

        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile

    def link(target_desc, objects,
             output_filename, output_dir=None, libraries=None,
             library_dirs=None, runtime_library_dirs=None,
             export_symbols=None, debug=0, extra_preargs=None,
             extra_postargs=None, build_temp=None, target_lang=None):

        print 'objects = ' + str(objects)

        for obj in objects[:]:
            if '_intel' in obj and '_gnu' in output_filename:
                objects.remove(obj)
            elif '_gnu' in obj and '_intel' in output_filename:
                objects.remove(obj)
            elif '_single' in obj and '_double' in output_filename:
                objects.remove(obj)
            elif '_double' in obj and '_single' in output_filename:
                objects.remove(obj)


        print 'objects = ' + str(objects)
        print 'target = ' + str(output_filename)

        if '-fopenmp' in extra_postargs:
            self.set_executable('linker_so', 'g++')
            self.set_executable('linker_exe', 'g++')

        elif '-openmp' in extra_postargs:
            self.set_executable('linker_so', 'icpc')
            self.set_executable('linker_exe', 'icpc')

        else:
            raise EnvironmentError('Unable to find a linker')
           
        super_link(target_desc, objects,
                   output_filename, output_dir, libraries,
                   library_dirs, runtime_library_dirs,
                   export_symbols, debug, extra_preargs,
                   extra_postargs, build_temp)

    self.link = link


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        build_ext.build_extensions(self)

print find_packages()

setup(name = pkg_name,
      version = '0.1',
      description = 'An open source method of characteristics code for ' + \
                  'solving the 2D neutron distribution in nuclear reactors]',
      author = 'Will Boyd',
      author_email = 'wboyd@mit.edu',
      url = 'https://github.com/mit-crpg/OpenMOC',

      # this is necessary so that the swigged python file gets picked up
#      py_modules=['openmoc'],
#      ext_package = 'openmoc',
      ext_modules = extensions,
      packages = find_packages(),
#      packages = ['openmoc'],
#      packages = ['openmoc', 'openmoc.intel.double', 'openmoc.intel.single', 
#                  'openmoc.gnu.double', 'openmoc.gnu.single'],
#      package_dir = {'openmoc.intel.double': 'openmoc/intel/double',
#                     'openmoc.intel.single': 'openmoc/intel/single',
#                     'openmoc.gnu.double': 'openmoc/gnu/double',
#                     'openmoc.gnu.single': 'openmoc/gnu/single',
#                     'openmoc': 'openmoc'},
#      packages = ['openmoc']

      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
