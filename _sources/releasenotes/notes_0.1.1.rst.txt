.. _notes_0.1.1:

===============================
Release Notes for OpenMOC 0.1.1
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

- `9ba06b2`_: Added exponential intrinsic evaluation option for GPUSolver.
- `0f2a8d8`_: Implemented a thread private vectorized solver class (VectorizedPrivateClass).
- `3ca79f5`_: Track segments now stored in a stride one array to improve cache coherency.
- `3baffee`_: OpenMP parallel regions use guided rather than dynamic scheduling.

.. _9ba06b2: https://github.com/mit-crpg/OpenMOC/commit/9ba06b2
.. _0f2a8d8: https://github.com/mit-crpg/OpenMOC/commit/0f2a8d8
.. _3ca79f5: https://github.com/mit-crpg/OpenMOC/commit/3ca79f5
.. _3baffee: https://github.com/mit-crpg/OpenMOC/commit/3baffee

---------
Bug Fixes
---------
