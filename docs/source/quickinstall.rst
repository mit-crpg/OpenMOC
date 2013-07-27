.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMOC on your computer. For more detailed instructions on configuring and installing OpenMOC, see :ref:`usersguide_install` in the User's Manual.

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMOC source code is hosted on GitHub_. To download and install OpenMOC, you need to install Git_ and the GNU C++ compiler_. In addition, you need to install Python_ 2.6 or 2.7 and the following Python packages on your machine: SWIG_, NumPy_, matplotlib_, h5py_ and setuptools_. These packages can easily be installed using a package manager for Linux and Mac OS (see :ref:`usersguide_install` for more details). The following command will install all required and optional dependencies on Ubuntu 12.04 or later::

    sudo apt-get install build-essential git swig python-numpy python-matplotlib python-h5py python-setuptools

If you have already installed each of these prerequisites, you can download and install OpenMOC by entering the following commands in the console::

    git clone git://github.com/mit-crpg/OpenMOC.git
    cd OpenMOC
    sudo python setup.py install

This will build a shared library accessible as a Python package named ``openmoc`` and install it (by default in /usr/local/lib/pythonX.X/dist-packages). The ``openmoc`` Python package can now be imported into any Python script as follows:

.. code-block:: python

    import openmoc


.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _Git: http://git-scm.com
.. _compiler: http://gcc.gnu.org/
.. _Python: http://www.python.org/
.. _SWIG: http://www.swig.org/
.. _NumPy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _h5py: http://www.h5py.org/
.. _setuptools: http://pythonhosted.org/setuptools/
