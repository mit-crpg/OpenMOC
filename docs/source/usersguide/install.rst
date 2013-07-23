.. _usersguide_install:

==============================
Installation and Configuration
==============================

-------------
Prerequisites
-------------

There are a few prerequisites which must be installed on your machine before you will be able to build and install OpenMOC. All of the prerequisites can easily be installed using a standard package manager (such as apt-get_ for Ubuntu, yum_ for Red Hat, and macports_ for Mac OS). The following instructions detail which packages are needed, with instructions on how to install them on Ubuntu.

.. admonition:: Required

    * GNU's C++ compiler_

      In order to compile OpenMOC, you will need to have a C++ compiler installed 
      on your machine. It is recommended that you build OpenMOC with g++ version 
      4.4 or later. 

      To install g++ on Debian Linux or a Debian derivative such as Ubuntu, 
      use the following command in the console::

	sudo apt-get install build-essential


    * Python_

      OpenMOC extensively uses Python for rapid data processing and visualization. 
      OpenMOC uses a mixture of fast, compiled C/C++/CUDA code with Python bindings 
      such that users are entirely insulated from having to write in C/C++ or CUDA, and
      can simply write a Python script with calls to OpenMOC.

      Nearly all common Linux distributions, as well as Mac OS, come with Python 
      pre-installed. Currently, OpenMOC is only tested with Python 2.6 and 2.7 and does
      not yet support Python 3 or greater.

      If you are not sure whether you have Python or not, you can check in the console
      by typing the following::

	python --version

      If you do not have Python, you can easily install version 2.7 in Ubuntu as follows::

	sudo apt-get install python2.7


    * Git_

      The OpenMOC development team uses the Git_ version control software for code
      development and collaboration. To download the code from our GitHub_ repository,
      you must first install git. You can do so on Ubuntu as follows::

	sudo apt-get install git


    * SWIG_

      The Simplified Wrapper Interface Generator (SWIG) is used to wrap the compiled 
      C/C++ and CUDA source code and generate Python bindings. SWIG can be installed on 
      Ubuntu as follows::
	
	sudo apt-get install swig


    * NumPy_

      NumPy is a commonly used, general purpose numerical array creation and 
      manipulation package for Python. OpenMOC uses NumPy to seamlessly transfer 
      numerical data between a Python script and the C/C++ back-end compiled code. To 
      install NumPy on Ubuntu, use the following command::

	sudo apt-get install python-numpy


.. admonition:: Optional

    * matplotlib_

      The matplotlib Python package is a general purpose utility for plotting and 
      visualizations. OpenMOC uses matplotlib in the ``openmoc.plotter`` module to 
      generate plots, such as flux and power distributions. Although matplotlib and
      the ``openmoc.plotter`` modules are not required for an OpenMOC simulation,
      visualizations can be a helpful tool in debugging and input model validation.

      To install matplotlib in Ubuntu, issue the following command::

	sudo apt-get install python-matplotlib


    * h5py_

      The h5py Python package contains tools to create, retrieve, and manipulate data 
      stored using the HDF5_ binary format. OpenMOC leverages the h5py functionality 
      in its ``openmoc.materialize`` module which includes routines to import/export 
      multi-group nuclear cross-section data to/from HDF5. In addition, the 
      ``openmoc.process`` module contains routines to import/export simulation and 
      performance data to/from HDF5 output files. Although h5py and the 
      ``openmoc.materialize`` and ``openmoc.process`` modules are not required to 
      run OpenMOC, they can facilitate fast access and storage of binary data.
      
      To install h5py on Ubuntu, issue the following command::
      
        sudo apt-get install python-h5py


.. _apt-get: http://www.apt-get.org/
.. _yum: http://yum.baseurl.org/
.. _macports: http://www.macports.org/
.. _compiler: http://gcc.gnu.org/
.. _Python: http://www.python.org/
.. _Git: http://git-scm.com
.. _SWIG: http://www.swig.org/
.. _NumPy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _h5py: http://www.h5py.org/
.. _HDF5: http://www.hdfgroup.org/HDF5/

--------------------
Obtaining the Source
--------------------

All OpenMOC source code is hosted on GitHub_. You can download the source code directly from GitHub or, if you have the Git_ version control software installed on your computer, you can use git to obtain the source code. The latter method has the benefit that it is easy to receive updates directly from the GitHub repository. GitHub has a good set of instructions_ for how to set up git to work with GitHub since this involves setting up ssh_ keys. With git installed and setup, the following command will download the full source code from the GitHub repository::

    git clone git://github.com/mit-crpg/OpenMOC.git

.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _Git: http://git-scm.com
.. _ssh: http://en.wikipedia.org/wiki/Secure_Shell
.. _instructions: http://help.github.com/set-up-git-redirect


--------------------
Building From Source
--------------------

Standard Build Configuration
----------------------------

To compile and install the code in a user local directory (recommended), simply run the following from the console::

  python setup.py install --user

To compile and install the code in the directory of all Python packages accessible to all users of your machine (not recommended), run the following command::

  python setup.py install

The code will now be accessible as a Python module from anywhere on your system.
The main OpenMOC Python package can be imported into any Python script as follows:

.. code-block:: python

    import openmoc


Custom Build Configuration
--------------------------

OpenMOC provides a number of user options to customize what and how OpenMOC source is compiled and built on your system. OpenMOC makes use of Python's distutils_ build configuration management module. 

To view a list of all of build commands supported by Python distutils, type the following in the console::
  
  python setup.py --help-commands

To install OpenMOC, we typically recommend using the :program:`install` command which builds and installs the code alongside other commonly referenced Python packagaes. The :program:`install` command includes its own set of options, some of which are defined by OpenMOC and some of which are defined by distutils_. To view a list of these options, type the following in the console::

  python setup.py install --help

The following section will provide an overview of the most useful and relevant build options for OpenMOC developers.

.. option:: --user

Installs OpenMOC in a user local directory (ie, `/home/username/.local/lib/pythonX.X/site-packages`) where it will only be accessible to your username. Installation without this option will instead install OpenMOC in the main Python directory accessible to all users of your machine (ie, `/usr/lib/pythonX.X/site-packages/`). This option is highly recommended for developers as it will prevent your Python packages from being polluted with code that has not yet been validated.


.. option:: --prefix=<path to install OpenMOC>

Installs OpenMOC to an explicitly defined directory. This options is generally not useful unless your directory is included in your :envvar:`PYTHONPATH` such that you can import ``openmoc`` into your Python scripts.


.. option:: --cc=<gcc,icpc,bgxlc>
	   
Sets the C++ compiler for the main ``openmoc`` module. Presently, GNU's gcc_, Intel's icpc_ and IBM's bgxlc_ are all configured if the path to the binary is pointed to by by the :envvar:`PATH` environment variable. The default setting is the :program:`gcc` compiler.


.. option:: --fp=<single,double>

Sets the floating point precision level for the main ``openmoc`` module. This sets the :envvar:`FP_PRECISION` macro in the source code by setting it as an environment variable at compile time. The default setting is :envvar:`single`.


.. option:: --with-cuda

Compiles the ``openmoc.cuda`` module using the :program:`nvcc` compiler. This module contains :cpp:class:`GPUSolver` class with MOC routines for execution on NVIDIA GPUs. The default build configuration does not include the ``openmoc.cuda`` module.


.. option:: --with-gcc

Compiles the ``openmoc.gnu.single`` and / or ``openmoc.gnu.double`` modules using GNU's :program:`gcc` C++ compiler. If one or both of :option:`--with-sp` and :option:`--with-sp` are also specified, the appropriate modules will be built and installed. If the floating point precision is not specified, ``openmoc.gnu.single`` will be built by default.


.. option:: --with-icpc

Compiles the ``openmoc.intel.single`` and / or ``openmoc.intel.double`` modules using Intel's :program:`icpc` C++ compiler. If one or both of :option:`--with-sp` and :option:`--with-sp` are also specified, the appropriate modules will be built and installed. If the floating point precision is not specified, ``openmoc.intel.single`` will be built by default.


.. option:: --with-bgxlc

Compiles the ``openmoc.bgxlc.single`` and / or ``openmoc.bgxlc.double`` modules using IBM's :program:`bgxlc` C++ compiler. If one or both of :option:`--with-sp` and :option:`--with-sp` are also specified, the appropriate modules will be built and installed. If the floating point precision is not specified, ``openmoc.bgxlc.single`` will be built by default.


.. option:: --with-sp

Specifies the :envvar:`single` floating point precision level to be used for ``openmoc.gnu.single``, ``openmoc.intel.single``, and / or ``openmoc.bgxlc.single`` modules. Compiles This option must be used in conjunction with the :option:`--with-gcc`, :option:`--with-icpc`, and / or :option:`--with-bgxlc` options.


.. option:: --with-dp

Specifies the :envvar:`double` floating point precision level to be used for ``openmoc.gnu.double``, ``openmoc.intel.double``, and / or ``openmoc.bgxlc.double`` modules. Compiles This option must be used in conjunction with the :option:`--with-gcc`, :option:`--with-icpc`, and / or :option:`--with-bgxlc` options.


.. option:: --debug-mode

Compiles with debugging symbols and information by including the :envvar:`-g` compile flag.


.. option:: --with-ccache

Compiles using ccache_ which uses a cache to speedup compilation of unchanged source files with the binaries from previous compilations. This flag is only relevant for developers needing to frequently recompile the source code. The ccache program must be installed for this flag to work.


.. option:: --with-papi

Compiles all :cpp:class:`Solver` derived classes with PAPI_ instrumentation for performance counter measurements.

.. _distutils: http://docs.python.org/2/library/distutils.html#module-distutils
.. _gcc: http://gcc.gnu.org/
.. _icpc: http://software.intel.com/en-us/intel-compilers
.. _bgxlc: http://www-03.ibm.com/software/products/us/en/ccompfami/
.. _ccache: http://ccache.samba.org
.. _NVIDIA: http://www.nvidia.com/content/global/global.php
.. _PAPI: http://icl.cs.utk.edu/papi/


-----------------------------
Installing on Ubuntu with PPA
-----------------------------

A binary package for Debian Linux derivatives, such as Ubuntu, is under development. Please check back at a later time for further updates.

