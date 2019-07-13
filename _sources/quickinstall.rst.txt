.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMOC on your computer. For more detailed instructions on configuring and installing OpenMOC, see :ref:`Installation <install>` in the User's Manual.

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMOC source code is hosted on GitHub_. To download and install OpenMOC, you need to install Git_ and the GNU C++ compiler_. In addition, you need to install Python_ version 3.1 or later and the following Python packages on your machine: SWIG_, NumPy_, matplotlib_, and h5py_. These packages can easily be installed using a package manager for Linux and Mac OS (see :ref:`Installation <install>` for more details). The following command will install all required and optional dependencies on Ubuntu 12.04 or later::

    sudo apt-get install build-essential git swig python-dev python-numpy python-matplotlib python-h5py

If you have already installed each of these prerequisites, you can download and install OpenMOC by entering the following commands in the console::

    git clone https://github.com/mit-crpg/OpenMOC.git
    cd OpenMOC
    python setup.py install --user

This will build a shared library accessible as a Python package named ``openmoc`` and install it (by default in /home/YourUserName/.local/lib/pythonX.X/dist-packages). The ``openmoc`` Python package can now be imported into any Python script as follows:

.. code-block:: python

    import openmoc

.. warning:: The :option:`--user` flag should be used verbatim and should **NOT** be replaced with your username.
.. warning:: Python 2.6 is no longer officially supported with OpenMOC, however the number of incompatibilities with Python 2 is currently very low.


.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _Git: http://git-scm.com
.. _compiler: http://gcc.gnu.org/
.. _Python: http://www.python.org/
.. _SWIG: http://www.swig.org/
.. _NumPy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _h5py: http://www.h5py.org/
.. _symbolic links: http://en.wikipedia.org/wiki/Symbolic_link
