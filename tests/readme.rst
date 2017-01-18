==================
OpenMOC Test Suite
==================

The purpose of this test suite is to ensure that that all user input options
can be used successfully without breaking the code. The test suite is based on
regression or integrated testing where different types of input files are
configured and the full OpenMOC code is executed. Results from simulations
are compared with expected results. The test suite is comprised of many tests
which reside in sub-directories in the tests directory.

The test suite is designed to integrate with cmake using ctest_. To run the
full test suite run:

.. code-block:: sh

    python run_tests.py

A subset of build configurations and/or tests can be run. To see how to use
the script run:

.. code-block:: sh

    python run_tests.py --help

As an example, say we want to run all tests with that have "adjoint" and
"pwr" in their name. Also, we wish to split the tests across 4 processors.
We can run:

.. code-block:: sh

    python run_tests.py -j 4 -R "adjoint|pwr"

Note that standard regular expression syntax is used for selecting build
configurations and tests. To print out a list of build configurations, we
can run:

.. code-block:: sh

    python run_tests.py -l

Note that the test suite requires h5py, matplotlib (>=1.5) and Pandas (>=0.13).

.. _ctest: http://www.cmake.org/cmake/help/v2.8.12/ctest.html
