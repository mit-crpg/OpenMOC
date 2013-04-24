import  os
from os.path import join as pjoin
from setuptools import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
import subprocess
import numpy

def find_in_path(name, path):
    "Find a file in a search path"
    #adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, or set $CUDAHOME')
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home':home, 'nvcc':nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in cudaconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))

    return cudaconfig
CUDA = locate_cuda()


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


openmoc = Extension('_openmoc',
                    sources=['openmoc/openmoc.i',
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
                             'openmoc/src/dev/DeviceMaterial.cu'],
                library_dirs=[CUDA['lib64']],
                swig_opts=['-c++', '-keyword'],
                libraries=['cudart'],
                runtime_library_dirs=[CUDA['lib64']],
                # this syntax is specific to this build system
                # we're only going to use certain compiler args with nvcc and not with gcc
                # the implementation of this trick is in customize_compiler() below
                extra_link_args=['-lstdc++', '-openmp', '-liomp5', '-lpthread',
                                 '-lirc', '-limf'],
                extra_compile_args={'gcc': ['-O3', '-openmp',  '-std=c++0x', '-fPIC'],
                                    'nvcc': ['-arch=sm_20', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                include_dirs = [numpy_include, CUDA['include'], 'openmoc/src/host', 'openmoc/src/dev'])


# check for swig
if find_in_path('swig', os.environ['PATH']):
    subprocess.check_call('swig -python -c++ -o openmoc/openmoc_wrap.cpp openmoc/openmoc.i', shell=True)
else:
    raise EnvironmentError('the swig executable was not found in your PATH')



def customize_compiler_for_nvcc(self):
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
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, openmoc, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', 'nvcc')
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            self.set_executable('compiler_so', 'icpc')
#            self.set_executable('linker_so', 'icpc')
            self.set_executable('linker_exe', 'icpc')
            postargs = extra_postargs['gcc']

        super(obj, src, openmoc, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)

setup(name='openmoc',
      # random metadata. there's more you can supploy
      author='Will Boyd',
      version='0.1',

      # this is necessary so that the swigged python file gets picked up
      py_modules=['openmoc'],
      package_dir={'': 'openmoc/src'},

      ext_modules = [openmoc],
      pagkages = ['openmoc'],

      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
