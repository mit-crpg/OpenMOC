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
``with_ccache``                    boolean                         False
``debug_mode``                     boolean                         False
``profile_mode``                   boolean                         False
``with_cuda``                      boolean                         False
``with_numpy``                     boolean                         True
``sources``                        dictionary of strings           C/C++/CUDA source files flags for each ``Extension`` object
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

  python setup.py install --user --fp=double

Similarly, the following command would build and install the ``openmoc`` module with default compiler (gcc), as well was the ``openmoc.cuda`` module, each with single precision and `debug symbols`_::
  
  python setup.py install --user --with-cuda --fp=single --debug-mode


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

The ``customize_compiler(...)`` method is used to override the ``_compile(...)`` method in distutils to allow for compilation with a variety of toolchains (*e.g.*, ``gcc``, ``icpc``, etc.). As presently implemented, the method chooses a compiler based on the macro definitions in the compile line (i.e, ``gcc`` for the macro definition :envvar:`-DGNU`). Likewise, the ``customize_linker(...)`` method is used to override the ``link(...)`` method in distutils to allow for linking with a variety of toolchains (*e.g.*, ``g++``, ``icpc``, etc.). The method chooses an executable for linking based on the target shared library name.

The ``custom_build_ext(...)`` class is used to override (subclass) the ``build_ext`` class in the ``distutils.command`` module. In particular, this class overrides the ``build_extension(...)`` method and uses it for the following:

- Inject the ``customize_compiler(...)`` and ``customize_linker(...)`` methods into the ``build_ext`` class
- Call SWIG to generate Python wrappers for the C/C++ source code.


Building and Installation
-------------------------

In the final step, the ``setup(...)`` method from the ``distutils.core`` module is called in the ``setup.py`` file. The ``setup(...)`` method receives the list of the ``Extension`` class objects and builds and installs each one as a shared library in the :file:`/home/<username>/.local/lib/python-x.x/site-packages/` directory. On a Unix-based machine, the shared library for the default ``openmoc`` module will be ``_openmoc.so``. The Python modules in OpenMOC (*e.g.*, ``openmoc.materialize``, ``openmoc.plotter``, etc.) will be installed in the :file:`/home/<username>/.local/lib/python-x.x/site-packages/` directory. 


--------------------
SWIG Interface Files
--------------------

OpenMOC uses the SWIG system (discussed in :ref:`Simplified Wrapper Interface Generator <swig>`) to generate Python bindings for classes and routines in the compiled C/C++ source code. In order for SWIG to work, the C/C++ header files **must contain all of the class and function prototypes.** Furthermore, the headers files must be exposed to SWIG through a `SWIG interface file`_ (see :ref:`SWIG Input <swig_input>`). OpenMOC includes the interface file :file:`/OpenMOC/openmoc/openmoc.i` for the main ``openmoc`` Python module, and an interface file :file:`/OpenMOC/openmoc/cuda/openmoc_cuda.i` for the ``openmoc.cuda`` module.

The :ref:`Add a C/C++ Source File <add_source_file>` section discusses how to add new C/C++ source files and expose them to SWIG through the interface files. The interface files are useful for a variety of auxiliary purposes as well, most notably the specifications to input and retrieve NumPy_ data from the compiled C/C++ shared library object(s) from Python (see :ref:`NumPy Typemaps <numpy_typemaps>`).


--------------------
Common Modifications
--------------------

The following sections cover some of the most common modifications and extensions that many developers need as they incorprate new features into OpenMOC, including how to add new C/C++ source files, Python modules, C/C++ extension modules and more.


.. _add_source_file:

Add a C/C++ Source File
-----------------------

There are three steps which must be taken to integrate a new source C/C++ file into the build system for OpenMOC. 
 
1. Include source header file (``.h``) in top of SWIG interface file (*e.g.* :file:`/OpenMOC/openmoc/openmoc.i`) using the following code syntax:

  .. code-block:: none
      
     #include "../src/MyFile.h"
    
2. Include source header file (``.h``) in bottom of SWIG interface file (*e.g.* :file:`/OpenMOC/openmoc/openmoc.i`) using the following code syntax:

  .. code-block:: none
		  
     %include ../src/MyFile.h

3. Append the source implementation file (``.c``, ``.cpp``, ``.cu``, etc.) to the ``sources`` attribute for the ``configuration`` class in the :file:`/OpenMOC/config.py` file (see :ref:`Configuration File <configuration_file>`).

.. note:: Changes to the C/C++ source files are not reflected until the OpenMOC has been reinstalled.


Add a Python Module
-------------------

OpenMOC includes several Python modules by default (i.e., ``openmoc.materialize``, ``openmoc.plotter``, etc.). These modules are Python files located in the :file:`/OpenMOC/openmoc` directory and are installed each time the C/C++ extension module(s) for OpenMOC are built and installed. For example, to create the ``openmoc.mymodule`` module, create the :file:`mymodule.py` file in the :file:`OpenMOC/openmoc` directory. You must then append the name of the module (*i.e.*, ``openmoc.mymodule``) to the ``packages`` list class attribute in the ``configuration`` class in the :file:`/OpenMOC/config.py` file.

.. note:: Changes to a Python module are not reflected until OpenMOC has been reinstalled.


Add a Compiler Flag
--------------------

In order to add a new compiler flag to OpenMOC, simply append it as a Python string to the ``compiler_flags`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``compiler_flags`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the compiler flag is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`-falign-functions` flag for ``gcc``, append '-falign-functions' to the list in ``compiler_flags`` corresponding to 'gcc'.


Add a Linker Flag
-----------------

In order to add a new linker flag to OpenMOC, simply append it as a Python string to the ``linker_flags`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``linker_flags`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the linker flag is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`-dynamiclib` flag for ``gcc``, append '-dynamiclib' to the list in ``linker_flags`` corresponding to 'gcc'.

Add a Library Directory
-----------------------

In order to add a new library directory to OpenMOC, simply append it as a Python string to the ``library_directories`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``library_directories`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the library directory is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`/usr/local/cuda/lib` include directory for ``nvcc``, append '/usr/local/cuda/lib' to the list in ``library_directories`` corresponding to 'gcc'.

You do not need to add a library directory if it is already in included in your :envvar:`LD_LIBRARY_PATH` environment variable. You can check if the directory is included with the following command from a Linux or Mac bash console:

.. code-block:: bash

    $ echo $LD_LIBRARY_PATH

To include the directory in :envvar:`LD_LIBRARY_PATH` instead of the :file:`config.py` file, use the following command from a Linux or Mac bash console:

.. code-block:: bash
      
    $ export LD_LIBRARY_PATH=/my/library/directory/here:$LD_LIBRARY_PATH


Add an Include Directory
------------------------

In order to add a new include directory to OpenMOC, simply append it as a Python string to the ``include_directories`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``include_directories`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the include directory is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`/usr/local/cuda/include` include directory for ``gcc``, append '/usr/local/cuda/include' to the list in ``include_directories`` corresponding to 'gcc'.


Link a Shared Library
---------------------

In order to link OpenMOC to a shared library, simply append the library name as a Python string to the ``shared_libraries`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``shared_libraries`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the shared library is applicable to and append the string to the list corresponding to that key. For example, to add the :option:`gomp` shared library for ``gcc``, append 'gomp' to the list in ``shared_libraries`` corresponding to 'gcc'.

.. warning:: Do **NOT** prepend "l" or "lib" to the shared library as is typical for most compilers. Python distutils will automatically do this for you.


Add a Macro Definition
----------------------

In order to add a C/C++ pre-processing macro option to OpenMOC, simply append the macro as a Python tuple to the ``macros`` attribute of the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The ``macros`` attribute is a Python dictionary (see :ref:`Configuration File <configuration_file>`) with keys for each compiler supported by the build system. Simply choose which compiler the macro is applicable to and append the tuple to the list corresponding to that key. For example, to add the :option:`METHOD=fast` macro for ``gcc``, append the ``('METHOD', 'fast')`` tuple to the list in ``macros`` corresponding to 'gcc'.


----------------------------
Add a C/C++ Extension Module
----------------------------

Many developers may write new C/C++ source code which performs some new physics or compute-intensive task. In some cases, the new code may be applicable for some users but less desirable for others. Alternatively, the code may only be applicable for certain types of simulations and less so for others. In these situations, it may be best to include the new C/C++ source code as a new extension module to OpenMOC. The following section discuses the steps which must be taken (in order) to incorporate a new extension module into OpenMOC's build system.


Create New SWIG Interface File
------------------------------

The first step to creating a new extension module is to create new SWIG interface file for the module. Interface files and some SWIG capabilities are discussed in more detail in :ref:`Simplified Wrapper Interface Generator <swig>`. In this section, it suffices to say that if you wish to create the ``openmoc.submodule`` extension module, you will need to create the :file:`openmoc_submodule.i` interface file akin to what is illustrated below:

.. code-block:: none

    %module openmoc_submodule

    %{ 

      #define SWIG_FILE_WITH_INIT

      /* Include all header files to wrap here */
      #include "first_source_file.h"
      ...
      #include "last_source_file.h"
    %}

    /* Include all header files to wrap here */
    %include first_source_file.h
    ...
    %include last_source_file.h


Add Source Files
----------------

Second, you need to create a new entry in the ``sources`` dictionary attribute for the ``configuration`` class in the :file:`/OpenMOC/config.py` file. The entry should contain a list of the C/C++ source files to compile, including the :ref:`SWIG wrap file <swig_input>`. An example of what might be appended to the ``configuration`` class is illustrated below

.. code-block:: python

    # Store list of C/C++ source files for the module
    sources['submodule'] = ['submodule_wrap.cpp',
		            'first_source_file.cpp',
			    ...
			    'last_source_file.cpp']

.. _add_build_options:

Add Build Options
-----------------

Next, you must create new entries in the ``configuration`` class attributes for compiler flags, linker flags, etc. in the :file:`/OpenMOC/config.py` file. Each of these build options is stored as a Python dictionary. Some of your build options may be identical to those for one or more of the main extension modules for OpenMOC. An example of what might be appended to the ``configuration`` class is illustrated below:

.. code-block:: python

    # Store build options for the module

    # Compiler flags
    compiler_flags['submodule'] = ['-first-option', ... , '-last-option']

    # Linker flags
    linker_flags['submodule'] = ['-first-option', ... , '-last-option']

    # Shared libraries (do not prepend "l" or "lib")
    shared_libraries['submodule'] = ['-firstlib', ..., '-lastlib']

    # Library directories (if not set in LD_LIBRARY_PATH)
    library_directories['submodule'] = ['/first/library/directory',
                                        ...
					'/second/library/directory']

    # Include directories
    include_directories['submodule'] = ['/first/include/directory',
                                        ...
					'/second/include/directory']

					
    # Define new macros as with tuples for each
    macros['submodule'] = [('FIRST_MACRO', None),
                           ...
			   ('LAST_MACRO', 42)]
    
Add a Build Option
------------------

In order to make your ``openmoc.submodule`` extension module an optional module for OpenMOC, you need to add a build option (or flag). The build options are defined and interpreted by the ``custom_install`` class in the :file:`/OpenMOC/setup.py` file. There are three primary steps which must be taken in order to add a build option.

1) First, append your option to the ``user_options`` list attribute in the ``custom_install`` class. The ``user_options`` list contains tuples of three elements each. The first element is the option string name, the second element is typically ``None``, and the third element is a brief description of the option to be printed for the ``--help`` option. An example of the creation of the :option:`--newopt` option which can be assigned a numerical or string value from the command line is given as follows:

   .. code-block:: python

      user_options = [...,
		      ('newopt=', None, 'A new option which takes a string or numerical value'),
		      ...]

2) Alternatively, if your option is a boolean option - for example, a binary switch to turn on/off the compilation of your ``openmoc.submodule`` extension module - you will need to define it in the ``boolean_options`` list attribute in the ``custom_install`` class:

   .. code-block:: python

       boolean_options = [...
                          'newopt',
			  ...]

3) Create a class attribute corresponding to your build option to the ``custom_install`` class in the :file:`/OpenMOC/setup.py` file. The class attribute should have the same name and capitalization as the command line option. This should be done in the ``initialize_options(...)`` class method. An example is given for the boolean option case as follows:

   .. code-block:: python

       def initialize_options(self):
       
         install.initialize_options(self)

	 ...
	 
	 Set a default value for the new build option
	 self.newopt = False

4) Inform the ``configuration`` class (given by the global ``config`` variable in ``setup.py``) to use the value passed in for the build option. This may entail creating a new class attribute for the ``configuration`` class in the :file:`/OpenMOC/config.py` file to account specifically for this option.

   .. code-block:: python

       def finalize_options(self):

         install.finalize_options(self)

	 ...

	 # Tell the configuration class what to do with the option
	 config.newopt = self.newopt

Once the build option is given to the ``configuration`` class, it is up to you to figure out what to do with it. The command line options are typically applied to the build configuration in the ``setup_extension_modules`` method for the ``configuration`` class in the :file:`/OpenMOC/config.py` file as presented in the next section.

	
Inject Compiler Options
-----------------------

Although you created a list of compiler flags for your module in :ref:`Add Build Options <add_build_options>`, they will not be used unless we expose them to the underlying distutils build system. We do this using the ``customize_compiler(...)`` routine in the :file:`/OpenMOC/setup.py` file. This routine overrides the ``_compile(...)`` routine in distutils and allows us to configure the compiler executable and compiler flags as we please. In particular, ``customize_compiler`` defines a new version of the ``_compile(...)`` method and *injects* it into distutils.

In order to ensure that distutils will use your compiler flags, you need to add a new conditional block to the ``_compile(...)`` method. This conditional should be able to determine when the ``_compile(...)`` method has received the build information for your extension module through its parameter set. One way to do this is to assign a specific macro to your build module which can be found in the ``pp_opts`` (pre-processing options) parameter to the ``_compile(...)`` method. for example, if you add a macro :envvar:`FIRST_MACRO` which is unique and only used for your extension module, then you can look for it in ``pp_opts`` and configure the compiler when needed. An example of how this might be done in the ``_compile(...)`` routine is shown as follows:

.. code-block:: python

    # If we find your macro, distutils is building your module
    if '-DFIRST_MACRO' in pp_opts:
	
      # Set the compiler executable to the compiler 
      # you want to use for your module
      self.set_executable('compiler_so', 'gcc')

      # Set the compiler flags
      postargs = config.compiler_flags['gcc']


Append Module to Installation Packages
--------------------------------------

Next, you must append the name of the module (*i.e.*, ``openmoc.submodule``) to the ``packages`` list class attribute in the ``configuration`` class in the :file:`/OpenMOC/config.py` file. This will ensure that your new extension module will not only be compiled, but also installed into the directory with the rest of the OpenMOC shared library extension module(s) and Python modules. An example is as follows: 

.. code-block:: python

    packages = ['openmoc', ..., 'openmoc.materialize', ..., 'openmoc.submodule']


Create SWIG Wrap File
---------------------

Your extension module C/C++ source files must be "wrapped" using SWIG to create a SWIG "wrap file." The distutils package will automatically do this for most Python distributions, but in some cases it is not done properly. To account for the inconsistencies across platforms, the OpenMOC build system manually calls SWIG for each extension module we wish to build at compile time. In particular, the ``custom_build_ext(...)`` routine in the :file:`/OpenMOC/setup.py` file is used to directly call SWIG to wrap the source code for each extension module. An example of what one might append to the ``custom_build_ext(...)`` routine to wrap the source for the ``openmoc.submodule`` extension module might be the following:

.. code-block:: python

    if 'submodule' in config.extensions:
      os.system('swig -python -c++ -keyword -o ' + \
		'openmoc/submodule/openmoc_submodule_wrap.cpp ' + \
                'openmoc/submodule/openmoc_submodule.i')


Create an Extension Object
--------------------------

The final step is to create an ``Extension`` object for your module in the ``setup_extension_modules(...)`` class method in the ``configuration`` class. The particular setup of your module is highly dependent on the functionality of your module and the build options added in the preceding steps. That said, there are several common issues to note when creating the ``Extension`` object.

First, the filename should begin with an underscore "_" since this is common practice for Python C/C++ extension modules. Second, the ``name`` parameter should be set to the filename of the shared library you wish to create. Typically, it is recommended that the shared library filename be constructed using the same words intended for the Python module, with "." replaced by "_". For example, to create the ``module.submodule.subsubmodule`` C/C++ extension module for Python, the shared library should be called :file:`_module_submodule_subsubmodule.so`. Finally, each of the parameters to the ``Extension`` object should be set using the lists of build options (*i.e.*, compiler flags, include libraries, etc.) configured in :ref:`Add Build Options <add_build_options>`.

An example of how one might instantiate the ``Extension`` object in the ``setup_extension_modules(...)`` routine is given below:

.. code-block:: python

    # Create Extension object and append to the list of objects
    # in the configuration class
    if config.submodule:

      self.extensions.append(
        Extension(name = '_openmoc_submodule', 
		  sources = self.sources['submodule'], 
		  library_dirs = self.library_directories['submodule'], 
		  libraries = self.shared_libraries['submodule'],
		  extra_link_args = self.linker_flags['submodule'], 
		  include_dirs = self.include_directories['submodule'],
		  define_macros = self.macros['submodule'][self.fp],
		  swig_opts = self.swig_flags,
		  export_symbols = ['init_openmoc']))

.. _python_package_tree:

Python Package Tree
-------------------

Finally, you will need to create a directory tree to represent your module within the OpenMOC `Python package`_. For example, for the ``openmoc.submodule`` extension module, you would need to create the :file:`/OpenMOC/openmoc/submodule` directory. In addition, you will need to create a Python :file:`__init__.py` file, which is required for Python to treat the directory as a Python package. For a C/C++ extension module, the :file:`__init__.py` file typically will import the shared library (*i.e.*, :file:`_openmoc_submodule.so`) as follows:

.. code-block:: python

    import _openmoc_submodule

    # Do any other initialization needed when
    # someone imports your module into Python
    ...

.. note:: You must create an :file:`__init__.py` file for each level of the Python package hierarchy. For example, for the ``openmoc.submodule.subsubmodule``, you will need to create the :file:`/OpenMOC/openmoc/submodule/__init__.py` and :file:`/OpenMOC/openmoc/submodule/subsubmodule/__init__.py` files.


Build the Extension Module
--------------------------

FINALLY, you should be prepared to compile and install your C/C++ extension module using OpenMOC's distutils-based build system!! The next step is to build OpenMOC as discussed in :ref:`Installation and Configuration <install>`, using the build option(s) which you defined for your module. For example, if you defined the build option :option:`--newopt` as a binary option for building your module, you might run the following in the console: 

.. code-block:: bash

    $ python setup.py install --user --newopt

This command should compile and install your shared library extension module in :file:`/home/<username>/.local/lib/pythonX.X/site-packages/` directory. In particular, you should be able to find the :file:`_openmoc_submodule.so` file in that directory, as well as the import directory tree for the module as presented in :ref:`Python Package Tree <python_package_tree>` (*i.e.*, :file:`/openmoc/submodule/` with the :file:`__init__.py` file).

If all was properly configured as described in the preceding steps, you should be able to import your extension module into Python with the following:

.. code-block:: python

    from openmoc.submodule import *

    # Do cool things with your extension module here
    ...

============================
Alternative C++ Build System
============================
Some developers might not wish to include the Python/SWIG build system due to the additional requirements such as SWIG_, Python_, Numpy_, matplotlib_, and h5py_. Since OpenMOC source code is entirely written in C++, it is possible to bypass the Python/SWIG build system and run OpenMOC using a compiled C++ input file. This is primarily useful for developers running OpenMOC on new architectures that might not support Python and SWIG. Performance analysis can also be easier without the Python/SWIG interface. 

.. note:: It is **highly** recommended that users utilize the regular Python/SWIG build system unless there is a specific reason for using the alternative C++ build system.

--------
Makefile
--------

The alternative C++ build system is available in the :file:`OpenMOC/profile/` 
directory. 
The build system depends on the source files found in the :file:`OpenMOC/src/` 
directory so that any changes to those files are noticed by both the regular 
Python/SWIG build system and the alternative C++ build system. 
The alternative C++ build system is based on the use of Make_. 
Since OpenMOC includes some C++11 features, developers should ensure 
the C++ compiler they wish to use includes full C++11 support. Developers using
the alternative C++ build system should be familiar with both C++ and Make_.

The alternative C++ build system compiler options are included in 
:file:`OpenMOC/profile/Makefile`. This Makefile includes several command line 
options that can be executed from the :file:`OpenMOC/profile/` directory:

  * **make** - Compiles the OpenMOC source files with the C++ input file
    indicated by the ``case`` variable in the Makefile. The default is 
    :file:`OpenMOC/profile/models/run_time_standard/run_time_standard.cpp`.
  * **make all** - Compiles the OpenMOC source files with all input
    files given in the Makefile indicated in the variable ``cases``.
  * **make run** - Runs the OpenMOC C++ input file indicated by the ``case``
    variable in the Makefile. Note: the case must be compiled before this
    command will run correctly.
  * **make clean** - Deletes all output files formed from compiling OpenMOC source 
    and input files described in the ``cases`` variable in the Makefile.

.. note:: Make_ does not handle all source dependencies, it is advised to delete the obj/ directory after changing any function's arguments or if using link time optimization

---------------------------------------
Building a multi-purpose C++ executable
---------------------------------------

The default case built with the Makefile is :file:`OpenMOC/profile/models/run_time_standard/run_time_standard.cpp`.
:file:`run_time_standard` is an OpenMOC executable that can read ".geo" geometry files, and takes the simulation parameters as command line input.
This can help users run large number of cases with the C++ build without recompiling.

-----------------------------------
Building individual C++ Input Files
-----------------------------------

It is advised to first compile the included example C++ inputs provided in the 
:file:`OpenMOC/profile/models/` directory. After ensuring the inputs compile and
verifying correct behavior runtime behavior, the developer should use the example
C++ inputs as a reference to write new C++ input files. Most OpenMOC commands 
should translate reasonably well between Python on C++ inputs, though some
aspects such as array declaration are much more difficult in C++. It is important
to note that ``Material``, ``Cell``, ``Universe``, and ``Lattice`` objects 
should be allocated on the heap rather than the stack to prevent segmentation 
faults when the program terminates. For instance, while

.. code-block:: cpp

    Cell basic_cell(0, "basic cell");
    Universe basic_universe(0, "basic universe");
    basic_universe.addCell(&basic_cell);
    Geometry geometry;
    geometry.setUniverses(1, 1, &basic_universe);

would seem like correct syntax, a segmentation fault would occur at the end of
the program execution when the geometry destructor is called since it deletes
all member universes. Likewise, universes delete all member cells. The correct
way of writing this code block to prevent this error would be to allocate
``basic_cell`` and ``basic_universe`` on the heap as:
 
.. code-block:: cpp

    Cell* basic_cell = new Cell(0, "basic cell");
    Universe* basic_universe = new Universe(0, "basic universe");
    basic_universe->addCell(basic_cell);
    Geometry geometry;
    geometry.setUniverses(1, 1, basic_universe);

----------------
Compiler Options
----------------

Once the C++ input file is completed, it should be added to the Makefile to be
compiled. To do this, add the location of the C++ input file to the ``cases``
variable in the Makefile. Note that the Makefile assumes input files are located
in a sub-directory of the :file:`OpenMOC/profile/models/` directory. Therefore, 
the Makefile adds :file:`models/` as a prefix to any input file locations.
In addition, several compiler options can be specified at the top of the 
Makefile, presented as variables detailed in :ref:`Table 2 <table_makefile_options>`.

.. _table_makefile_options:

=========================  ===============================================  ==============================   =========================
Variable                   Description                                      Allowed Values                   Default Value
=========================  ===============================================  ==============================   =========================
``COMPILER``               Specifies the compiler                           gnu, intel, clang, bluegene      mpicc
``OPENMP``                 Flag to turn on OpenMP parallelism               yes, no                          yes
``OPTIMIZE``               Flag to turn on compiler optimizations           yes, no                          yes
``DEBUG``                  Flag to turn on vector reports and debug mode    yes, no                          no
``PROFILE``                Creates a :file:`gmon.out` file for profiling    yes, no                          no
``INFO``                   Display information on compiler optimization     yes, no                          no
``PRECISION``              Specifies the floating point precision           single, double                   single
``CMFD_PRECISION``         Specifies the floating point precision in CMFD   single, double                   single
=========================  ===============================================  ==============================   =========================

**Table 2**: Makefile compiler options for the alternative C++ build system.


After the input is compiled, an executable will be created in the same directory
as the input file. The executable can be called directly based on its location
or it can be run using the **make run** command described previously with a 
modification of the ``case`` variable in the Makefile.


Geometry files (.geo)
---------------------

For complex geometries, it may be preferable to load the geometry in Python first and dump it to a ".geo" file, which can then be read by a C++ executable.

.. code-block:: python

    # Dump geometry to file
    geometry.dumpToFile("path_to_file")

.. code-block:: C++

    // Load geometry from .geo file
    geometry.loadFromFile("path_to_file");

.. _distutils: http://docs.python.org/2/library/distutils.html
.. _make: http://www.gnu.org/software/make/
.. _cmake: http://www.cmake.org/
.. _C/C++ Extensions: http://docs.python.org/2/extending/building.html
.. _shared library: http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
.. _debug symbols: http://en.wikipedia.org/wiki/Debug_symbol
.. _SWIG interface file: http://www.swig.org/Doc2.0/SWIGDocumentation.html#Introduction_nn6
.. _NumPy: http://www.numpy.org
.. _Python package: http://docs.python.org/2/tutorial/modules.html#packages
.. _SWIG: http://www.swig.org/
.. _NumPy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _h5py: http://www.h5py.org/
.. _Make: http://www.gnu.org/software/make/
