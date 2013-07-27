.. _usersguide_processing:

=================================
Data Processing and Visualization
=================================

This section is intended to explain in detail the recommended procedures for
carrying out common tasks with OpenMC. While several utilities of varying
complexity are provided to help automate the process, in many cases it will be
extremely beneficial to do some coding in Python to quickly obtain results.  In
these cases, and for many of the provided utilities, it is necessary for your
Python installation to contain:

* [1]_ `Numpy <http://www.numpy.org/>`_
* [1]_ `Scipy <http://www.scipy.org/>`_
* [2]_ `h5py <http://code.google.com/p/h5py/>`_
* [3]_ `Matplotlib <http://matplotlib.org/>`_

Most of these are easily obtainable in Ubuntu through the package manager, or
are easily installed with distutils.

.. [1] Required for array-based data manipulation and processing
.. [2] Required only if reading HDF5 materials data files
.. [3] Optional for plotting utilities

Most of these are easily obtainable in Ubuntu through the package manager, or
are easily installed with distutils.


----------------------
Geometry Visualization
----------------------


------------------
Flux Visualization
------------------
