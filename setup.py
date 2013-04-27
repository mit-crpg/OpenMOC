import os
from os.path import join as pjoin
from setuptools import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
import subprocess
import numpy

# Adapted from Robert McGibbon's CUDA distutils setup provide in open source 
# form here: https://github.com/rmcgibbo/npcuda-example


cpp_compiler = 'icpc'
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


sources = ['openmoc/src/host/Cell.cpp',
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

libraries = []
library_dirs = []
runtime_library_dirs = []
extra_link_args = []
extra_compile_args = []
include_dirs = []
extra_link_args = []
extra_compile_args = {}
include_dirs = [numpy_include, 'openmoc/src/host', 'openmoc/src/dev']
swig_opts = ['-c++', '-keyword']

if cpp_compiler == 'gcc':
    extra_link_args.extend(['-lstdc++', '-lgomp', '-fopenmp', '-shared'])
    extra_compile_args['gcc'] = ['-c', '-O3', '-fopenmp', '-std=c++0x', '-fPIC']

elif cpp_compiler == 'icpc':
    library_dirs.extend([ICPC['lib64']])
    runtime_library_dirs.extend([ICPC['lib64']])
    extra_link_args.extend(['-lstdc++', '-openmp', '-liomp5', '-lpthread', 
                            '-lirc', '-limf', '-lrt', '-lpython2.7', '-shared'])
    extra_compile_args['icpc'] =['-c', '-O3', '-openmp', '-std=c++0x', '-fpic',
                               '-xhost', '-openmp-report', '-vec-report']
    include_dirs.extend([ICPC['include']])


if use_cuda:
    sources.extend(['openmoc/src/dev/DeviceMaterial.cu',
                    'openmoc/src/dev/DeviceTrack.cu',
                    'openmoc/src/dev/DeviceFlatSourceRegion.cu',
                    'openmoc/openmoc-cuda.i'])
    libraries.extend(['cudart'])
    library_dirs.extend([CUDA['lib64']])
    runtime_library_dirs.extend([CUDA['lib64']])
    extra_compile_args['nvcc'] = ['-c', '-arch=sm_20', '--ptxas-options=-v', 
                                  '--compiler-options', '-fpic']
    include_dirs.extend([CUDA['include']])

else:
    sources.extend(['openmoc/openmoc.i'])


openmoc = Extension(name = '_openmoc', 
                    sources = sources, 
                    library_dirs = library_dirs, 
                    libraries = libraries,
                    runtime_library_dirs = runtime_library_dirs,
#                    extra_compile_args = extra_compile_args,
                    extra_link_args = extra_link_args, 
                    include_dirs = include_dirs,
                    swig_opts = swig_opts)


# check for swig
if find_in_path('swig', os.environ['PATH']):
    subprocess.check_call('swig -python -c++ -keyword -o ' + \
                          'openmoc/openmoc_wrap.cpp openmoc/openmoc.i', 
                          shell=True)
else:
    raise EnvironmentError('The swig executable was not found in your PATH')


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
#    self.linker_so = ['icpc', '-lstdc++', '-openmp', '-liomp5', \
#                      '-lpthreads', '-lirc', '-limf']
    super_compile = self._compile
    super_link = self.link

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, openmoc, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cpp' and cpp_compiler == 'gcc':
            self.set_executable('compiler_so', 'gcc')
            postargs = extra_compile_args['gcc']

        elif os.path.splitext(src)[1] == '.cpp' and cpp_compiler == 'icpc':
            self.set_executable('compiler_so', 'icpc')
            postargs = extra_compile_args['icpc']

        elif os.path.splitext(src)[1] == '.cu':
            self.set_executable('compiler_so', 'nvcc')
            postargs = extra_compile_args['nvcc']

            for item in cc_args:
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

        if cpp_compiler == 'gcc':
            self.set_executable('linker_so', 'g++')
            self.set_executable('linker_exe', 'g++')

        elif cpp_compiler == 'icpc':
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


setup(name = 'openmoc',
      version = '0.1',
      description = 'An open source method of characteristics code for ' + \
                  'solving the 2D neutron distribution in nuclear reactors]',
      author = 'Will Boyd',
      author_email = 'wboyd@mit.edu',
      url = 'https://github.com/mit-crpg/OpenMOC',

      # this is necessary so that the swigged python file gets picked up
      py_modules=['openmoc'],
      ext_modules = [openmoc],
      packages = ['openmoc'],

      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
