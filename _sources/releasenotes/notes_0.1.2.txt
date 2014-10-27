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

- `bd9d4c11eb`_: CMFD coarse group solver
- `ea25db0454`_: Python 3 compatibility
- `75731d8820`_: Single first-time installation

.. _bd9d4c11eb: https://github.com/mit-crpg/OpenMOC/commit/bd9d4c11eb
.. _ea25db0454: https://github.com/mit-crpg/OpenMOC/commit/ea25db0454
.. _75731d8820: https://github.com/mit-crpg/OpenMOC/commit/75731d8820


---------
Bug Fixes
---------

- `216ce92b2a`_: Minor correction to C5G7 cross-sections
- `ac6bbe7674`_: Fix to find coordinates within the ``Geometry`` in ``openmoc.plotter.plot_fluxes(...)`` routine
- `ea25db0454`_: Fix to allow sources to be stored in ``openmoc.process.store_simulation_state(...)`` routine
- `c9bd55788e`_: Segfault for coarse track spacings in ``Solver::convergeSource(...)`` routine
- `25fe960cd5`_: Fixed ``openmoc`` import statement in ``openmoc.process`` submodule

.. _216ce92b2a: https://github.com/mit-crpg/OpenMOC/commit/216ce92b2a
.. _ea25db0454: https://github.com/mit-crpg/OpenMOC/commit/ea25db0454
.. _ac6bbe7674: https://github.com/mit-crpg/OpenMOC/commit/ac6bbe7674
.. _c9bd55788e: https://github.com/mit-crpg/OpenMOC/commit/c9bd55788e
.. _25fe960cd5: https://github.com/mit-crpg/OpenMOC/commit/25fe960cd5
