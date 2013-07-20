.. _quickinstall:

===================
Quick Install Guide
===================

This quick install guide outlines the basic steps needed to install OpenMOC on
your computer. For more detailed instructions on configuring and installing
OpenMOC, see :ref:`usersguide_install` in the User's Manual.

..
.. --------------------------------
.. Installing on Ubuntu through PPA
.. --------------------------------
..
.. For users with Ubuntu 11.10 or later, a binary package for OpenMC is available through a `Personal Package Archive`_ (PPA) and can be installed through the `APT package manager`_. Simply enter the following commands into the terminal:

.. code-block:: sh

..    sudo apt-add-repository ppa:paulromano/staging
..    sudo apt-get update
..    sudo apt-get install openmc

.. _Personal Package Archive: https://launchpad.net/~paulromano/+archive/staging
.. _APT package manager: https://help.ubuntu.com/community/AptGet/Howto

-------------------------------------------
Installing from Source on Linux or Mac OS X
-------------------------------------------

All OpenMOC source code is hosted on GitHub_. To download and install OpenMOC, you need to install git_ and GNU's C++ compiler, g++_. In addition, you need to install Python 2.6 or 2.7 and the following Python packages on your machine: swig_, setuptools_, h5py_, numpy_, and matplotlib_. These packages can easily be installed using a package manager for Linux and Mac OS (see :ref:`usersguide_install` for more details).

You can download and install OpenMOC by entering the following commands in a terminal:

.. code-block:: sh

    git clone git://github.com/mit-crpg/OpenMOC.git
    cd OpenMOC
    sudo python setup.py install

This will build a shared library accessible as a Python package named ``openmoc`` and install it (by default in /usr/local/lib/python2.7/dist-packages). The ``openmoc`` Python package can now be imported into any Python script as follows:

.. code-block:: python

    import openmoc

.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _git: http://git-scm.com
.. _g++: http://gcc.gnu.org/
.. _swig: http://www.swig.org/
.. _h5py: http://www.h5py.org/
.. _numpy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _setuptools: https://pythonhosted.org/setuptools/index.html
