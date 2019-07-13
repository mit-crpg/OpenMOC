.. _notes_0.1.4:

===============================
Release Notes for OpenMOC 0.1.4
===============================

-------------------
System Requirements
-------------------

There are no special requirements for running the OpenMOC code. As of this
release, OpenMOC has been tested on a variety of Linux distributions as well as
Mac OS X. However, it has not been tested yet on any releases of Microsoft
Windows. Memory requirements will vary depending on the size of the problem at
hand (the number of flat source regions and the level of track discretization).

------------
New Features
------------

- `751ba0`_: Major reformulation of constructive solid geometry models used in OpenMOC
- `691afb`_: Inclusion of (n,2n), (n,3n), etc. reaction with nu-scattering matrix
- `41ac75`_: Compatiblity module for the soon-to-be-open-sourced OpenCG package
- `e69bbf`_: Directory parameter for store_simulation_state(...) routine

.. _751ba0 : https://github.com/mit-crpg/OpenMOC/commit/751ba09225cdec74168c4d43fb34848c2668ad97
.. _691afb : https://github.com/mit-crpg/OpenMOC/commit/691afb2fa71a861e21c2e2a8635e741e35845109
.. _41ac75 : https://github.com/mit-crpg/OpenMOC/commit/41ac756973f6c3770aa4820d2e4e68e02da7bc9e
.. _e69bbf : https://github.com/mit-crpg/OpenMOC/commit/e69bbfd0d59a6f2c124762f30f1f7d557b509c18

---------
Bug Fixes
---------

- `0936be`_: Revised import statements in Python modules to be compatible with Python 3
- `189733`_: Bug fixes to VectorizedSolver class caused by backwards compatiblity issues from v0.1.3
- `732017`_ : Bug fix for openmoc.compatible.casmo module's k-infinity parser

.. _0936be : https://github.com/mit-crpg/OpenMOC/commit/0936bec595e88a423e86e5b4e84accd32b11e647
.. _189733 : https://github.com/mit-crpg/OpenMOC/commit/1897dd76ce0821cc2038477ab5b814de205bb602
.. _732017 : https://github.com/mit-crpg/OpenMOC/commit/732017226be65857b6d838393e5a000d11a79583
