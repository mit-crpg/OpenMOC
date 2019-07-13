=======
OpenMOC
=======


.. image:: https://api.travis-ci.org/mit-crpg/OpenMOC.svg?branch=develop
    :target: https://travis-ci.org/mit-crpg/OpenMOC
.. image:: https://coveralls.io/repos/github/mit-crpg/OpenMOC/badge.svg?branch=develop
    :target: https://coveralls.io/github/mit-crpg/OpenMOC?branch=develop
.. image:: https://img.shields.io/badge/powered%20by-OpenMOC-blue.svg
    :target: https://mit-crpg.github.io/OpenMOC/
.. image:: https://img.shields.io/badge/pypi-v0.4.0-orange.svg
    :target: https://pypi.python.org/pypi/openmoc/0.4.0
.. image:: https://img.shields.io/badge/license-MIT%20License-brightgreen.svg    
    :target: https://mit-crpg.github.io/OpenMOC/license.html
.. image:: https://img.shields.io/badge/anucene-Elsevier-lightgray.svg
    :target: http://www.sciencedirect.com/science/article/pii/S0306454913006634

Welcome to the OpenMOC repository! OpenMOC is a simulation tool for
solving for the flux, power distribution, and multiplication factor
within a nuclear reactor. The code employs the deterministic method
of characteristics, with support for both fixed source and eigenvalue
calculations. The OpenMOC project aims to provide a simple-to-use
Python package bound to a back-end of source code written in C/C++
and CUDA. It includes support for constructive solid geometry and 2D
ray tracing for fully heterogeneous multi-group calculations.
Development of OpenMOC began at MIT in 2012 and is spearheaded by
several graduate students in the
`Nuclear Science & Engineering Department`_.

Complete documentation on OpenMOC is hosted at
https://mit-crpg.github.io/OpenMOC/. If you would like to
contribute to the OpenMOC project, please `contact`_ the
development team.

For a guided example, see a demonstration `IPython Notebook`_.

------------
Installation
------------

Detailed `installation instructions`_ can be found in the
User's Guide.

---------------
Troubleshooting
---------------

Join the OpenMOC `users group`_ to ask questions and discuss
methods and simulation workflows.

--------------
Citing OpenMOC
--------------

Please cite OpenMOC in your publications if it helps your research:

.. code-block:: latex

    @article{openmoc2014,
      author = {Boyd, William and Shaner, Samuel and Li, Lulu and Forget, Benoit and Smith, Kord},
      journal = {Annals of Nuclear Energy},
      title = {The OpenMOC Method of Characteristics Neutral Particle Transport Code},
      volume = {68},
      pages = {43--52},
      year = {2014}
    }

-------
License
-------

OpenMOC is approved for distribution under the MIT/X license_.

.. _installation instructions: https://mit-crpg.github.io/OpenMOC/usersguide/install.html
.. _license: https://mit-crpg.github.io/OpenMOC/license.html
.. _Nuclear Science & Engineering Department: http://web.mit.edu/nse/
.. _IPython Notebook: http://nbviewer.ipython.org/gist/anonymous/abbce6824bceda49a615
.. _contact: https://mit-crpg.github.io/OpenMOC/developers.html
.. _users group: https://groups.google.com/forum/#!forum/openmoc-users
