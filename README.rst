=================================================================
OpenMOC Method of Characteristics Neutral Particle Transport Code
=================================================================

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

------------
Installation
------------

Detailed `installation instructions`_ can be found in the 
User's Guide.

---------------
Troubleshooting
---------------

If you run into problems installing or running OpenMOC, 
please post your issue to the Google Group `users forum`_. 

--------
Citation
--------

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
.. _users forum: https://groups.google.com/forum/#!forum/openmoc-users
