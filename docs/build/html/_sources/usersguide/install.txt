.. _usersguide_install:

==============================
Installation and Configuration
==============================

--------------------
Building from Source
--------------------

Prerequisites
-------------


In addition, you need to install Python 2.6 or 2.7 and the following Python packages on your machine: swig_, setuptools_, h5py_, numpy_, and matplotlib_. These packages can easily be installed using a package manager for Linux and Mac OS (see :ref:`usersguide_install` for more details).

These packages can be installed from the console using a package manager for Ubuntu, Fedora and Mac OS as follows:


Ubuntu:

.. code-block:: sh

    sudo apt-get install swig python-setuptools python-numpy python-matplotlib python-h5py

Fedora / Red Hat:

.. code-block:: sh

    sudo yum install swig python-setuptools numpy python-matplotlib
    sudo easy_install python-h5py

Mac OS:

.. code-block:: sh

    sudo port install python26
    sudo port install swig-python py26-numpy py26-matplotlib
    sudo easy_install python-setuptools h5py


.. admonition:: Required

    * GNU's g++ C++ Compiler

      In order to compile OpenMOC, you will need to have a C++ compiler
      installed on your machine. It is recommended that you build OpenMOC
      with g++ version 4.4.6 or later. 

      To install g++ on Debian Linux or a Debian derivative such as Ubuntu, 
      use the following command in the console::

	sudo apt-get install build-essential

    * git_ Version Control Software::

	sudo apt-get install git



    * Simplified Wrapper Interface Generator (swig_)

      SWIG is used to wrap the compiled C/C++ and CUDA source code and
      generate Python bindings. SWIG can be installed on Debian derivatives
      as follows::
	
	sudo apt-get install swig

    * Python's setuptools_ Package

      Installation and configuration management is handled by the 
      setuptools_ Python package. To install setuptools_, use the 
      following command::

	sudo apt-get install python-setuptools

    * Python's numpy_ Package

      The numpy_ Python package is a commonly used, general purpose 
      numerical array creation and manipulation. OpenMOC uses numpy_
      to seamlessly transfer numerical data from a Python script to
      the C/C++ back-end compiled code. To install numpy_, use the
      following command::

	sudo apt-get install python-numpy


.. admonition:: Optional

    * Python's h5py_ Package

      The h5py_ Python package contains tools to create, retrieve, and 
      manipulate data stored using the HDF5_ binary format. OpenMOC 
      leverages the h5py_ functionality in its ``openmoc.materialize``
      module which includes routines to import/export multi-group
      materials nuclear cross-section data to/from HDF5_. In addition,
      the ``openmoc.process`` module contains routines to import/export
      simulation and performance data to/from HDF5_ output files. 
      Although h5py_ and the ``openmoc.materialize`` and ``openmoc.process``
      modules are not required to run OpenMOC, they can facilitate fast
      storage and access of binary data.
      
      To install h5py_ on a Debian derivative, issue the following command::
      
        sudo apt-get install python-h5py

    * Python's matplotlib_ Package

      The matplotlib_ Python package is a general purpose utility for 
      plotting and visualizations. OpenMOC uses matplotlib_ in the 
      ``openmoc.plotter`` module to generate plots, such as flux
      and power distributions. Although the matplotlib_ and the 
      ``openmoc.process`` module is not required for an OpenMOC simulation,
      visualizations can be a helpful tool in debugging and input model
      validation.

      To install matplotlib_ from Debian Linux, issue the following command
      in the terminal::

	sudo apt-get install python-matplotlib


.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _git: http://git-scm.com
.. _g++: http://gcc.gnu.org/
.. _swig: http://www.swig.org/
.. _h5py: http://www.h5py.org/
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _numpy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _setuptools: https://pythonhosted.org/setuptools/index.html

Obtaining the Source
--------------------

All OpenMOC source code is hosted on GitHub_. You can download the source code
directly from GitHub or, if you have the git_ version control software installed
on your computer, you can use git to obtain the source code. The latter method
has the benefit that it is easy to receive updates directly from the GitHub
repository. GitHub has a good set of `instructions
<http://help.github.com/set-up-git-redirect>`_ for how to set up git to work
with GitHub since this involves setting up ssh_ keys. With git installed and
setup, the following command will download the full source code from the GitHub
repository::

    git clone git://github.com/mit-crpg/OpenMOC.git

.. _git: http://git-scm.com
.. _ssh: http://en.wikipedia.org/wiki/Secure_Shell


Building From Source
--------------------

To compile and install the code in a user local directory (recommended), simply run the following from the console::

  python setup.py install --user

To compile and install the code in the directory of all Python packages accessible to all users of your machine (not recommended), run the following command::

  python setup.py install

The code will now be accessible as a Python module from anywhere on your system.
The main OpenMOC Python package can be imported into any Python script as follows:

.. code-block:: python

    import openmoc

