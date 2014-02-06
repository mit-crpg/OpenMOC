.. _install:

==============================
Installation and Configuration
==============================

-------------
Prerequisites
-------------

* :ref:`Installing Prerequisites on Ubuntu <ubuntu_prerequisites>`
* :ref:`Installing Prerequisites on Mac OS X <mac_prerequisites>`


--------------------
Obtaining the Source
--------------------

All OpenMOC source code is hosted on GitHub_. You can download the source code directly from GitHub or, if you have the Git_ version control software installed on your computer, you can use git to obtain the source code. The latter method has the benefit that it is easy to receive updates directly from the GitHub repository. GitHub has a good set of instructions_ for how to set up git to work with GitHub since this involves setting up ssh_ keys. With git installed and setup, the following command will download the full source code from the GitHub repository::

    git clone https://github.com/mit-crpg/OpenMOC.git

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

.. warning:: The :option:`--user` flag should be used verbatim and should **NOT** be replaced with your username.

To compile and install the code in the directory of all Python packages accessible to all users of your machine (not recommended), run the following command::

  python setup.py install

The code will now be accessible as a Python module from anywhere on your system.
The main OpenMOC Python package can be imported into any Python script as follows:

.. code-block:: python

    import openmoc

.. warning:: Some Python distributions require that the code be built and installed *twice* in order to setup all of the `symbolic links`_ correctly. This is only the case the first time OpenMOC is installed on a machine. Hence, code developers who subsequently make modifications to the source code only need to build and install the code once to have access to the modified binary.


Custom Build Configuration
--------------------------

OpenMOC provides a number of user options to customize what and how OpenMOC source is compiled and built on your system. OpenMOC makes use of Python's distutils_ build configuration management module. 

To view a list of all of build commands supported by Python distutils, type the following in the console::
  
  python setup.py --help-commands

To install OpenMOC, we typically recommend using the :program:`install` command which builds and installs the code alongside other commonly referenced Python packagaes. The :program:`install` command includes its own set of options, some of which are defined by OpenMOC and some of which are defined by distutils_. To view a list of these options, type the following in the console::

  python setup.py install --help


.. _build_configuration_options:

Build Configuration Options
---------------------------

This section section will provides an overview of the most useful and relevant build options for OpenMOC developers.

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

Compiles the ``openmoc.bgq.single`` and / or ``openmoc.bgq.double`` modules using IBM's :program:`bgxlc` C++ compiler. If one or both of :option:`--with-sp` and :option:`--with-sp` are also specified, the appropriate modules will be built and installed. If the floating point precision is not specified, ``openmoc.bgq.single`` will be built by default.


.. option:: --with-sp

Specifies the :envvar:`single` floating point precision level to be used for ``openmoc.gnu.single``, ``openmoc.intel.single``, and / or ``openmoc.bgq.single`` modules. This option must be used in conjunction with the :option:`--with-gcc`, :option:`--with-icpc`, and / or :option:`--with-bgxlc` options.


.. option:: --with-dp

Specifies the :envvar:`double` floating point precision level to be used for ``openmoc.gnu.double``, ``openmoc.intel.double``, and / or ``openmoc.bgq.double`` modules. This option must be used in conjunction with the :option:`--with-gcc`, :option:`--with-icpc`, and / or :option:`--with-bgxlc` options.


.. option:: --debug-mode

Compiles with debugging symbols and information by including the :envvar:`-g` compile flag.


.. option:: --no-numpy

Compiles OpenMOC without embedding the NumPy C API. This is severely limiting for integrating both OpenMOC source convergence calculations and data analysis into Python, but may be necessary on some machines such as IBM's BlueGene_ where NumPy is not a standard package.


.. option:: --with-papi

Compiles all :cpp:class:`Solver` derived classes with PAPI_ instrumentation for performance counter measurements.


.. option:: --with-ccache

Compiles using ccache_ which uses a cache to speedup compilation of unchanged source files with the binaries from previous compilations. This flag is only relevant for developers needing to frequently recompile the source code. The ccache p]rogram must be installed for this flag to work. The following console command will install ccache on Ubuntu::

    sudo apt-get install ccache


.. _distutils: http://docs.python.org/2/library/distutils.html#module-distutils
.. _gcc: http://gcc.gnu.org/
.. _icpc: http://software.intel.com/en-us/intel-compilers
.. _bgxlc: http://www-03.ibm.com/software/products/us/en/ccompfami/
.. _ccache: http://ccache.samba.org
.. _NVIDIA: http://www.nvidia.com/content/global/global.php
.. _PAPI: http://icl.cs.utk.edu/papi/
.. _symbolic links: http://en.wikipedia.org/wiki/Symbolic_link
.. _BlueGene: http://www-03.ibm.com/systems/technicalcomputing/solutions/bluegene/

-----------------------------
Installing on Ubuntu with PPA
-----------------------------

A binary package for Debian Linux derivatives, such as Ubuntu, is under development. Please check back at a later time for further updates.

