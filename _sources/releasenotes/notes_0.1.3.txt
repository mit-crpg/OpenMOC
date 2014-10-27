.. _notes_0.1.3:

===============================
Release Notes for OpenMOC 0.1.3
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

- `bee1e60397`_: Reconstructed segmentation routines to allow for tracking across overlapping cells.
- `a5ccd1324f`_: Added P0 boundary flux update to CMFD solver.

.. _a5ccd1324f : https://github.com/mit-crpg/OpenMOC/commit/a5ccd1324f
.. _bee1e60397 : https://github.com/mit-crpg/OpenMOC/commit/bee1e60397

---------
Bug Fixes
---------

- `a5ccd1324f`_: Fixed procedural error in computing diffusion coefficients in CMFD.
- `bee1e60397`_: Fixed error in computing the CMFD optically thick correction factor.
- `bee1e60397`_: TrackGenerator geometry string was incorrect resulting in TrackGenerator often not identifying a valid track file when one actually existed. The error in generating the geometry string has been corrected.
- `bee1e60397`_: Fixed inconsistency in XPlane and YPlane constants.

.. _a5ccd1324f : https://github.com/mit-crpg/OpenMOC/commit/a5ccd1324f
.. _bee1e60397 : https://github.com/mit-crpg/OpenMOC/commit/bee1e60397
