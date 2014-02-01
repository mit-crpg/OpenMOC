.. _build_sytem:

===================
Python Build System
===================

OpenMOC is built and installed using the Python distutils_ package, which is a standard package provided with all Python distributions. Some developers will find this to be an adjustment since this stands in contrast to the legacy make_ and cmake_ build systems used by many scientific simulation codes. The following sections will describe how a developer can interact with the build system for configuration management.


.. _configuration_file:

------------------
Configuration File
------------------

The primary options for the build and installation of OpenMOC are defined in the :file:`/OpenMOC/config.py` file. In particular, :file:`config.py` creates the `C/C++ Extensions`_ using distutils in the ``configuration`` class. The ``Extension`` class in the ``distutils.extension`` module defines all of the source code and compiler and linker options to create a `shared library`_ which can be imported into Python as module. The class attributes for the ``configuration`` class encapsulate the build options specified by the user at compile time (see :ref:`Build Configuration Options <build_configuration_options>`) and are itemized in :ref:`Table 1 <table_configuration_attributes>`.

.. _table_configuration_attributes:

=================================  ==============================  ===========================================  
Class Attribute                    Type                            Default Value 
=================================  ==============================  ===========================================
``cc``                             string                          'gcc'                                        
``fp``                             string                          'single'
``cpp_compilers``                  list of strings                 []
``fp_precision``                   list of string                  []
``with_ccache``                    boolean                         False
``debug_mode``                     boolean                         False
``with_cuda``                      boolean                         False
``with_numpy``                     boolean                         True
``sources``                        dictionary of strings           C/C++/CUDA source files flags for eac  ``Extension`` object
``compiler_flags``                 dictionary of strings           Flags for each ``Extension`` object
``linker_flags``                   dictionary of strings           Flags for each ``Extension`` object
``shared_libraries``               dictionary of strings           Libraries for each ``Extension`` object
``library_directories``            dictionary of strings           Libraries for each ``Extension`` object
``include_directories``            dictionary of strings           Includes for each ``Extension`` object
``macros``                         dictionary of strings           Macros for each ``Extension`` object
``extensions``                     list of ``Extension`` objects   
=================================  ==============================  ===========================================

**Table 1**: Attributes for the ``configuration`` class.

The class attributes each have default values corresponding to the configuration needed to build and install the ``openmoc`` module. The :file:`OpenMOC/setup.py` file instantiates a ``configuration`` class object and redefines the class attributes based on the command line options given at runtime. This process is discussed in greater detail in the :ref:`Setup File <setup_file>` section.

Once the build options have been defined for the ``configuration`` class, the :file:`setup.py` calls the ``setup_extension_modules(...)`` class method which instantiates one or more ``Extension`` class objects for each Python module the user wishes to build and install. 

For example, the following command would build and install the ``openmoc`` module using the default compiler (gcc) and with double precision::

  python setup.py install --user --with-dp

Users may build and install more than one OpenMOC module at once with the appropriate command line options. For example, the following command would build and install the ``openmoc`` module with default compiler (gcc), as well as the ``openmoc.intel.double`` module, each with double precision::

  python setup.py install --user --with-intel --with-dp

Similarly, the following command would build and install the ``openmoc`` module with default compiler (gcc), as well was the ``openmoc.cuda.single`` module, each with single precision and `debug symbols`_::
  
  python setup.py install --user --with-cuda --with-sp --debug-mode


.. _setup_file:


----------
Setup File
----------

The :file:`/OpenMOC/setup.py` file is used to execute the build and installation process for OpenMOC. The methodology and implementation to handle command line options, custom compilers and linkers, and the build and installation phase.


Command Line Options
--------------------

The :file:`setup.py` file first defines the ``custom_install`` class to override (subclass) the ``install`` class in the ``distutils.command.install`` module. In particular, the ``initialize_options(...)`` class method is used to initialize the default build options for the ``configuration`` class (see :ref:`Configuration File <configuration_file>`). The ``finalize_options(...)`` routine extracts the command line options, initializes the corresponding ``configuration`` class attributes, and calls the ``configuration.setup_extension_modules()`` class method to create the distutils ``Extension`` object.


Custom Compilers and Linkers
----------------------------

The ``customize_compiler(...)`` method is used to override the ``_compile(...)`` method in distutils to allow for compilation with a variety of toolchains (e.g., ``gcc``, ``icpc``, etc.). As presently implemented, the method chooses a compiler based on the macro definitions in the compile line (i.e, ``gcc`` for the macro definition :envvar:`-DGNU`). Likewise, the ``customize_linker(...)`` method is used to override the ``link(...)`` method in distutils to allow for linking with a variety of toolchains (e.g., ``g++``, ``icpc``, etc.). The method chooses an executable for linking based on the target shared library name.

The ``custom_build_ext(...)`` class is used to override (subclass) the ``build_ext`` class in the ``distutils.command`` module. In particular, this class overrides the ``build_extension(...)`` method and uses it for the following:

- Inject the ``customize_compiler(...)`` and ``customize_linker(...)`` methods into the ``build_ext`` class
- Call SWIG to generate Python wrappers for the C/C++ source code.


Building and Installation
-------------------------

In the final step, the ``setup(...)`` method from the ``distutils.core`` module is called in the ``setup.py`` file. The ``setup(...)`` method receives the list of the ``Extension`` class objects and builds and installs each one as a shared library in the :file:`/home/<username>/.local/lib/python-x.x/site-packages/` directory. On a Unix-based machine, the shared library for the default ``openmoc`` module will be ``_openmoc.so``. The Python modules in OpenMOC (e.g., ``openmoc.materialize``, ``openmoc.plotter``, etc.) will be installed in the :file:`/home/<username>/.local/lib/python-x.x/site-packages/` directory. 


.. _swig_interface_file:

--------------------
SWIG Interface Files
--------------------

OpenMOC uses the SWIG system (discussed in :ref:`Simplified Wrapper Interface Generator <swig>`) to generate Python bindings for classes and routines in the compiled C/C++ source code. In order for SWIG to work, the C/C++ header files **must contain all of the class and function prototypes.** Furthermore, the headers files must be exposed to SWIG through a `SWIG interface file`_. The interface files are located in the :file:`/OpenMOC/openmoc/...` directory and use a ``.i`` extension. There are different interface files for the different C/C++ extension modules which may be built for Python (e.g. with different compilers). :ref:`Table 2 <table_openmoc_swig_files>` tabulates all of the interface files in OpenMOC, the Python module that would be built, and the shell command that would be used to build the module.

.. _table_openmoc_swig_files:

============================================================  ========================  ================================================================  
File                                                          Python Module             Shell Build Command                                            
============================================================  ========================  ================================================================
:file:`/OpenMOC/openmoc/openmoc.i`                            ``openmoc``               :command:`python setup.py install --user`                       
:file:`/OpenMOC/openmoc/gnu/single/openmoc_gnu_single.i`      ``openmoc.gnu.single``    :command:`python setup.py install --user --with-gcc -with-sp`  
:file:`/OpenMOC/openmoc/gnu/double/openmoc_gnu_double.i`      ``openmoc.gnu.double``    :command:`python setup.py install --user --with-gcc --with-dp` 
:file:`/OpenMOC/openmoc/intel/single/openmoc_intel_single.i`  ``openmoc.intel.single``  :command:`python setup.py install --user --with-icpc --with-sp`
:file:`/OpenMOC/openmoc/intel/double/openmoc_intel_double.i`  ``openmoc.intel.double``  :command:`python setup.py install --user --with-icpc --with-dp` 
:file:`/OpenMOC/openmoc/bgq/single/openmoc_bgq_single.i`      ``openmoc.bgq.single``    :command:`python setup.py install --user --with-bgxlc --with-sp`
:file:`/OpenMOC/openmoc/bgq/double/openmoc_bgq_double.i`      ``openmoc.bgq.double``    :command:`python setup.py install --user --with-bgxlc --with-dp`
:file:`/OpenMOC/openmoc/cuda/single/openmoc_cuda.i`           ``openmoc.cuda``          :command:`python setup.py install --user --with-cuda`           
:file:`/OpenMOC/penmoc/cuda/single/openmoc_cuda_single.i`     ``openmoc.cuda.single``   :command:`python setup.py install --user --with-cuda --with-sp`
:file:`/OpenMOC/openmoc/cuda/double/openmoc_cuda_double.i`    ``openmoc.cuda.double``   :command:`python setup.py install --user --with-cuda --with-dp` 
============================================================  ========================  ================================================================ 

**Table 2**: SWIG interface files for OpenMOC modules. 

The :ref:`Add a C/C++ Source File <add_source_file>` section discusses how to add new C/C++ source files and expose them to SWIG through the interface files. The interface files are useful for a variety of auxiliary purposes as well, most notably the specifications to input and retrieve NumPy_ data from the compiled C/C++ shared library object(s) from Python (see :ref:`NumPy Typemaps <numpy_typemaps>`).


--------------------
Common Modifications
--------------------


.. _add_source_file:

Add a C/C++ Source File
-----------------------

There are three steps which must be taken to integrate a new source C/C++ file into the build system for OpenMOC. 
 
1. Include source header file (``.h``) in top of SWIG interface file (e.g. :file:`/OpenMOC/openmoc/openmoc.i`) using the following code syntax:

  .. code-block:: bash
      
     #include "../src/MyFile.h"
    
2. Include source header file (``.h``) in bottom of SWIG interface file (e.g. :file:`/OpenMOC/openmoc/openmoc.i`) using the following code syntax:

  .. code-block:: bash
		  
     %include ../src/MyFile.h

3. Append the source implementation file (``.c``, ``.cpp``, ``.cu``, etc.) to the ``sources`` attribute for the ``configuration`` class in the :file:`/OpenMOC/config.py` file (see :ref:`Configuration File <configuration_file>`).

.. note:: Changes to the C/C++ source files are not reflected until the OpenMOC has been reinstalled.


Add a Python Module
-------------------

OpenMOC includes several Python modules by default (i.e., ``openmoc.materialize``, ``openmoc.plotter``, etc.). These modules are Python files located in the :file:`/OpenMOC/openmoc` directory and are installed each time the C/C++ extension module(s) for OpenMOC are built and installed. For example, to create the :file:`mymodule.py` module, simply locate the file in the :file:`OpenMOC/openmoc` directory.

.. note:: Changes to a Python module are not reflected until the OpenMOC has been reinstalled.


Add a Compiler Flag
--------------------

In order to add a new compiler flag to OpenMOC, simply append it as a Python string to the ``compiler_flags`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``compiler_flags`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the compiler flag is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`-falign-functions` flag for ``gcc``, append '-falign-functions' to the list in ``compiler_flags`` corresponding to 'gcc'.


Add a Linker Flag
-----------------

In order to add a new linker flag to OpenMOC, simply append it as a Python string to the ``linker_flags`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``linker_flags`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the linker flag is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`-dynamiclib` flag for ``gcc``, append '-dynamiclib' to the list in ``linker_flags`` corresponding to 'gcc'.


Add an Include Directory
------------------------

In order to add a new include directory to OpenMOC, simply append it as a Python string to the ``include_directories`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``include_directories`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the include directory is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`/usr/local/cuda/include` include directory for ``gcc``, append '/usr/local/cuda/include' to the list in ``include_directories`` corresponding to 'gcc'.


Link a Shared Library
---------------------

In order to link OpenMOC to a shared library, simply append the library name as a Python string to the ``shared_libraries`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``shared_libraries`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the shared library is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`gomp` shared library for ``gcc``, append 'gomp' to the list in ``shared_libraries`` corresponding to 'gcc'.


Add a Macro Definition
----------------------

In order to add a C/C++ pre-processing macro option to OpenMOC, simply append the macro as a Python tuple to the ``macros`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``macros`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the macro is applicable to and append the tuple to the list corresponding to that key. For example, to add the :option:`METHOD=fast` macro for ``gcc``, append the ``('METHOD', 'fast')`` tuple to the list in ``macros`` corresponding to 'gcc'.



Add a C/C++ Extension Module
----------------------------



.. _distutils: http://docs.python.org/2/library/distutils.html
.. _make: http://www.gnu.org/software/make/
.. _cmake: http://www.cmake.org/
.. _C/C++ Extensions: http://docs.python.org/2/extending/building.html
.. _shared library: http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
.. _debug symbols: http://en.wikipedia.org/wiki/Debug_symbol
.. _SWIG interface file: http://www.swig.org/Doc2.0/SWIGDocumentation.html#Introduction_nn6
.. _NumPy: http://www.numpy.org
