.. _notes_0.1.2:

===============================
Release Notes for OpenMOC 0.1.2
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

- Reconstructed segmentation routines to allow for tracking across overlapping cells.
- Added P0 boundary flux update to CMFD solver.
- Increased parallel performance of CPUSolver.
- CMFD corner splitting now splits the corner currents to both neighboring surfaces with a weight of 0.5 instead of giving the current to only one surface.

---------
Bug Fixes
---------

- Fixed procedural error in computing diffusion coefficients in CMFD.
- Fixed error in computing the CMFD optically thick correction factor.
- TrackGenerator geometry string was incorrect resulting in TrackGenerator often not identifying a valid track file when one actually existed. The error in generating the geometry string has been corrected.
- Fixed inconsistency in XPlane and YPlane constants.
