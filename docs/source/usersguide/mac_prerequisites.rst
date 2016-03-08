.. _mac_prerequisites:

====================================
Installing Prerequisites on Mac OS X
====================================

There are a few prerequisites which must be installed on your machine before you will be able to build and install OpenMOC. All of the prerequisites can easily be installed using a standard package manager, such as MacPorts_ for Mac OS X. The following instructions detail which packages are needed, with instructions on how to install them using MacPorts.

.. note:: It is important to install the following prerequisite packages in the order given - in particular, Python must be installed before NumPy and Matplotlib.

.. admonition:: Required

    * GNU's C++ compiler_

      In order to compile OpenMOC, you will need to have a C++ compiler installed on your machine. It is recommended that you build OpenMOC with g++ version 4.4 or later. To install the latest stable version of ``g++`` using MacPorts (v4.8 at the time of this writing), use the following command in the console::

	sudo port install gcc48

      Finally, you must set this version of ``g++`` as the default for Macports::

	sudo port select --set gcc mp-gcc48


    * Python_

      OpenMOC extensively uses Python for rapid data processing and visualization. OpenMOC uses a mixture of fast, compiled C/C++/CUDA code with Python bindingssuch that users are entirely insulated from having to write in C/C++ or CUDA, and can simply write a Python script with calls to OpenMOC.

      Currently, OpenMOC has been tested with Python versions 2.6, 2.7, 3.1 and 3.2. Even if you already have Python installed, it is recommended that you install it again using MacPorts to ensure that it is properly configured with ``g++``. You can easily install version 2.7 using MacPorts as follows::

	sudo port install python27

      Finally, you must set this version of ``python`` as the default for Macports::

	sudo port select --set python python27


    * Git_

      The OpenMOC development team uses the Git_ version control software for code development and collaboration. To download the code from our GitHub_ repository, you must first install git. You can do so using MacPorts as follows::

	sudo port install git


    * SWIG_

      The Simplified Wrapper Interface Generator (SWIG) is used to wrap the compiled C/C++ and CUDA source code and generate Python bindings. SWIG can be installed using MacPorts::
	
	sudo port install swig swig-python


    * NumPy_

      NumPy is a commonly used, general purpose numerical array creation and manipulation package for Python. OpenMOC uses NumPy to seamlessly transfer numerical data between Python and the C/C++ back-end compiled code. To install NumPy using MacPorts, use the following command::

	sudo port install py27-numpy

.. admonition:: Optional

    * matplotlib_

      The matplotlib Python package is a general purpose utility for plotting and visualizations. OpenMOC uses matplotlib in the ``openmoc.plotter`` module to generate plots, such as flux and power distributions. Although matplotlib and the ``openmoc.plotter`` modules are not required for an OpenMOC simulation, visualizations can be a helpful tool in debugging and input model validation. To install matplotlib using MacPorts issue the following command::

	sudo port install py27-matplotlib


    * h5py_

      The h5py Python package contains tools to create, retrieve, and manipulate data stored using the HDF5_ binary format. OpenMOC leverages the h5py functionality in its ``openmoc.materialize`` module which includes routines to import/export multi-group nuclear cross-section data to/from HDF5. In addition, the ``openmoc.process`` module contains routines to import/export simulation and performance data to/from HDF5 output files. Although h5py and the ``openmoc.materialize`` and ``openmoc.process`` modules are not required to run OpenMOC, they can facilitate fast access and storage of binary data.
      
      To install h5py using MacPorts, issue the following command::
      
        sudo port install py27-h5py


.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _MacPorts: http://www.macports.org/
.. _compiler: http://gcc.gnu.org/
.. _Python: http://www.python.org/
.. _Git: http://git-scm.com
.. _SWIG: http://www.swig.org/
.. _NumPy: http://www.numpy.org/
.. _BlueGene: http://www-03.ibm.com/systems/technicalcomputing/solutions/bluegene/
.. _matplotlib: http://matplotlib.org/
.. _h5py: http://www.h5py.org/
.. _HDF5: http://www.hdfgroup.org/HDF5/
