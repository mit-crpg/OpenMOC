#! /usr/bin/env python

# System imports
import os
from distutils.sysconfig import get_config_vars
from distutils.core import *
from distutils      import sysconfig, dir_util

# Third-party modules - we depend on numpy for everything
import numpy


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

(opt,) = get_config_vars('OPT')
os.environ['OPT'] = ' '.join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)

ld = sysconfig._config_vars['LDSHARED']
sysconfig._config_vars['LDSHARED'] = ld.replace('gcc','icpc')

# range extension module
openmoc = Extension('_openmoc',
                   include_dirs=[numpy_include],
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
                            'openmoc/src/host/Universe.cpp'],
                   extra_compile_args=['-O3', '-fopenmp', '-std=c++0x'],
                   extra_link_args=['-lstdc++', '-lgomp', '-fopenmp'],
                   language='c++',
                   swig_opts=['-c++', '-keyword'],
)

#intel = ['-openmp', '-lirc', '-limf', '-liomp5', '-lpthread']

# NumyTypemapTests setup
setup(  name        = 'OpenMOC',
        description = 'An open source method of characteristics code for solving the 2D neutron distribution in nuclear reactor cores',
        author      = 'Will Boyd',
        author_email = 'wboyd@mit.edu',
        url = 'https://github.com/mit-crpg/OpenMOC',
        version     = '0.1',
        ext_modules = [openmoc],
        packages = ['openmoc'],
)
