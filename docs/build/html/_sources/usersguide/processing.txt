.. _usersguide_processing:

=================================
Data Processing and Visualization
=================================

This section is intended to explain in detail the recommended procedures for carrying out common tasks with OpenMOC. While several utilities of varying complexity are provided to help automate the process, in many cases it will be extremely beneficial to do some coding in Python to quickly obtain results. In these cases, and for many of the provided utilities, it is necessary for your Python installation to contain:

* `Numpy <http://www.numpy.org/>`_ - Required for array-based data manipulation and processing
* `h5py <http://www.h5py.org/>`_ - Required only if reading/writing HDF5 data files
* `Matplotlib <http://matplotlib.org/>`_ - Optional for plotting utilities


Most of these are easily obtainable in Ubuntu through the package manager.

-------------------------
Exporting Simulation Data
-------------------------

OpenMOC's ``openmoc.process`` module provides the ``storeSimulationState(...)`` routine to export simulation data to binary output files. The only required parameter for the routine is a ``Solver`` object. Optional parameters may be used to indicate whether to store the data in HDF5 or as a Python pickle_ file (default), store the fluxes, sources, pin powers and more. All of the supported parameters are listed in :ref:`Table 1 <table_store_simulation_state>`, and the output variables stored in the binary file are tabulated in :ref:`Table 2 <table_output_variables>`.

.. _table_store_simulation_state:

==============  ==================  ===================  ==========  ====================================
Parameter       Type                Default              Optional    Note
==============  ==================  ===================  ==========  ====================================
``solver``      ``Solver`` object   None                 No
``fluxes``      boolean             False                Yes         Whether to store the FSR fluxes
``sources``     boolean             False                Yes         Whether to store the FSR sources
``pin_powers``  boolean             False                Yes         Whether to store the pin powers
``use_hdf5``    boolean             False (pickle file)  Yes         Whether to use HDF5 
``filename``    string              'simulation-state'   Yes         The filename for storage
``append``      boolean             True                 Yes         Append to a file or create a new one
``note``        string              None                 Yes         Any comment on the simulation
==============  ==================  ===================  ==========  ====================================

**Table 1**: Parameters for the ``openmoc.proces.storeSimulationState(...)`` routine.

.. _table_output_variables:

=========================  ==============  =========================================
Output Variable            Type            Note
=========================  ==============  =========================================
solver type                string          'CPUSolver', 'ThreadPrivateSolver', etc.
# FSRs                     integer         
# materials                integer
# energy groups            integer
# tracks                   integer
# segments                 integer
track spacing [cm]         float
# azimuthal angles         integer
# polar angles             integer
# iterations               integer
source residual threshold  float
exponential                string          'exp intrinsic' or 'linear interpolation'
floating point             string          'double' or 'single'
time [sec]                 float           Total time to converge the source
keff                       float
note                       string          If requested by user
# threads                  integer         For solvers on multi-core CPUs
# threads per block        integer         For solvers on GPUs
FSR scalar fluxes          float array     If requested by user
FSR sources                float array     If requested by user
fission rates              float array(s)  If requested by user
=========================  ==============  =========================================

**Table 2**: Output variables in a binary file created by the ``openmoc.proces.storeSimulationState(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.process as proc

    # Setup and run simulation
    ...

    # Export the simulation data to an output file
    proc.storeSimulationState(solver, use_hdf5=True)



--------------------------
Retrieving Simulation Data
--------------------------

Exporting simulation data is only useful if there is a straightforward means to retrieve it for data processing at a later time. OpenMOC's ``restoreSimulationStates(...)`` routine in the ``openmoc.process`` module can be used for this purpose. This routine takes a binary data file created by the ``storeSimulationState(...)`` routine, parses the file and catalogues the data in a Python dictionary_, which it returns to the user. The parameters accepted by the routine are described in :ref:`Table 3 <table_restore_simulation_states>`, while the dictionary keys are identical to the output variables given in :ref:`Table 2 <table_output_variables>`.

.. _table_restore_simulation_states:

==============  =======  ======================  ========
Parameter       Type     Default                 Optional
==============  =======  ======================  ========
``filename``    string   'simulation-state.pkl'  Yes
``directory``   string   'simulation-states'     Yes
==============  =======  ======================  ========

**Table 3**: Parameters for the ``openmoc.process.restoreSimulationStates(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.process as proc

    # Retrieve the simulation state(s) stored in the 'states.h5' file
    # and returns the data in a Python dictionary
    simulation_state = proc.restoreSimulationState(filename='states.h5')


--------------------
Computing Pin Powers
--------------------

In some cases, a user may wish to only compute and export the pin powers for a simulation. In this case, the ``computeFSRPinPowers(...)`` routine in the ``openmoc.process`` module  may be used. The routine takes in a ``Solver`` subclass (e.g., ``ThreadPrivateSolver``, ``GPUSolver``, etc.) and computes the fission rate for each universe in the geometry by summing up the fission rates in each cell in the universe. In most cases, a universe is replicated in many places throughout the geometry. To account for this, the routine will separately compute the fission rates for each unique placement of that universe in the geometry. By default, the pin powers will be exported to a Python pickle_ file, but may alternatively be exported to an HDF5 binary file. :ref:`Table 4 <table_pin_powers>` describes the parameters accepted by the routine.

.. _table_pin_powers:

============  ==================  ========  =========
Parameter     Type                Default   Optional
============  ==================  ========  =========
``solver``    ``Solver`` object   None      No
``use_hdf5``  boolean             False     Yes
============  ==================  ========  =========

**Table 4**: Parameters for the ``openmoc.process.computeFSRPinPowers(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.process as proc

    # Setup and run simulation
    ...

    # Compute and export the pin powers
    proc.computeFSRPinPowers(solver, use_hdf5=True)

.. note:: The pin powers are computed for each nested universe level in the hierarchical geometry model.
.. note:: The pin powers are NOT normalized in any way - this is left to the user's discretion during data processing.


----------------------
Geometry Visualization
----------------------


Plotting Tracks
---------------

To plot the tracks crossing the geometry, use the ``plotTracks(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 5 <table_plot_tracks>`.

.. _table_plot_tracks:

===================  =========================  =========  =========  =========================
Parameter            Type                       Default    Optional   Note
===================  =========================  =========  =========  =========================
``track_generator``  ``TrackGenerator`` object  None       No         The tracks of interest
===================  =========================  =========  =========  =========================

**Table 5**: Parameters for the ``openmoc.plotter.plotTracks(...)`` routine.

The code snippet below illustrates the use of this routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geometry and generate tracks
    ...

    plot.plotTracks(geometry)

A depiction of the tracks for the :file:`/OpenMOC/sample-input/large-lattice.py` example input file with 4 azimuthal angles and 0.1 cm track spacing is illustrated in :ref:`Figure 1 <figure_tracks>`.

.. _figure_tracks:

.. figure:: ../../img/tracks.png
   :align: center
   :figclass: align-center
   :width: 400px

   **Figure 1**: The tracks crossing a a 4 :math:`\times` 4 lattice.

.. note:: The runtime required by the plotting routine scales with the number of tracks, which is proportional to the number of azimuthal angles and inversely proportional the track spacing.


Plotting Segments
-----------------

To plot the segments crossing the geometry color-coded by flat source region, use the ``plotSegments(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 6 <table_plot_segments>`.

.. _table_plot_segments:

===================  =========================  =========  =========  =========================
Parameter            Type                       Default    Optional   Note
===================  =========================  =========  =========  =========================
``track_generator``  ``TrackGenerator`` object  None       No         The tracks of interest
===================  =========================  =========  =========  =========================

**Table 6**: Parameters for the ``openmoc.plotter.plotSegments(...)`` routine.

The code snippet below illustrates the use of this routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geometry and generate tracks
    ...

    plot.plotSegments(geometry)

A depiction of the segments for the :file:`/OpenMOC/sample-input/large-lattice.py` example input file with 4 azimuthal angles and 0.1 cm track spacing is illustrated in :ref:`Figure 2 <figure_segments>`.

.. _figure_segments:

.. figure:: ../../img/segments.png
   :align: center
   :figclass: align-center
   :width: 400px

   **Figure 2**: The segments crossing a a 4 :math:`\times` 4 lattice.

.. warning:: This routine will require a long time for large geometries or fine track discretization. In addition, the Matplotlib consumes a substantial amount of memory to plot the segments and may throw a `segmentation fault`_ for large geometries.
.. note:: The runtime required by the plotting routine scales with the number of segments, which is proportional to the number of flat source regions and number of azimuthal angles and inversely proportional the track spacing.


Plotting by Material
--------------------

To plot the geometry color-coded by the material ID's throughout the geometry, use the ``plotMaterials(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 7 <table_plot_materials>`.

.. _table_plot_materials:

============  ===================  =========  =========  =========================
Parameter     Type                 Default    Optional   Note
============  ===================  =========  =========  =========================
``geometry``  ``Geometry`` object  None       No         The geometry of interest
``gridsize``  integer              250        Yes        The pixel resolution
============  ===================  =========  =========  =========================

**Table 7**: Parameters for the ``openmoc.plotter.plotMaterials(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geomery
    ...

    # Plot a 500 x 500 pixel image of the materials
    plot.plotMaterials(geometry, gridsize=500)

A depiction of the materials for the :file:`/OpenMOC/sample-input/large-lattice.py` example input file is illustrated in :ref:`Figure 3 <figure_materials>`.

.. _figure_materials:

.. figure:: ../../img/materials.png
   :align: center
   :figclass: align-center
   :width: 400px

   **Figure 3**: A 4 :math:`\times` 4 lattice color-coded by material.

.. note:: The runtime required by the plotting routine scales with the number of pixels in the image (the square of the ``gridsize`` parameter).
.. note:: The routine randomly selects a colormap at runtime. As a result, the colors in the figure will vary from run to run.


Plotting by Cell
----------------
To plot the geometry color-coded by the cell ID's throughout the geometry, use the ``plotCells(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 8 <table_plot_cells>`.

.. _table_plot_cells:

============  ===================  =========  =========  =========================
Parameter     Type                 Default    Optional   Note
============  ===================  =========  =========  =========================
``geometry``  ``Geometry`` object  None       No         The geometry of interest
``gridsize``  integer              250        Yes        The pixel resolution
============  ===================  =========  =========  =========================

**Table 8**: Parameters for the ``openmoc.plotter.plotCells(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geomery
    ...

    # Plot a 500 x 500 pixel image of the cells
    plot.plotCells(geometry, gridsize=500)

A depiction of the cells for the :file:`/OpenMOC/sample-input/large-lattice.py` example input file is illustrated in :ref:`Figure 4 <figure_cells>`.

.. _figure_cells:

.. figure:: ../../img/cells.png
   :align: center
   :figclass: align-center
   :width: 400px

   **Figure 4**: A 4 :math:`\times` 4 lattice color-coded by cell.

.. note:: The runtime required by the plotting routine scales with the number of pixels in the image (the square of the ``gridsize`` parameter).
.. note:: The routine randomly selects a colormap at runtime. As a result, the colors in the figure will vary from run to run.


Plotting by FSR
---------------

To plot the geometry color-coded by the flat source region ID's throughout the geometry, use the ``plotFlatSourceRegions(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 9 <table_plot_fsrs>`.

.. _table_plot_fsrs:

============  ===================  =========  =========  =========================
Parameter     Type                 Default    Optional   Note
============  ===================  =========  =========  =========================
``geometry``  ``Geometry`` object  None       No         The geometry of interest
``gridsize``  integer              250        Yes        The pixel resolution
============  ===================  =========  =========  =========================

**Table 9**: Parameters for the ``openmoc.plotter.plotFlatSourceRegions(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geomery
    ...

    # Plot a 500 x 500 pixel image of the flat source regions
    plot.plotFlatSourceRegions(geometry, gridsize=500)

A depiction of the flat source regions for the :file:`/OpenMOC/sample-input/large-lattice.py` example input file is illustrated in :ref:`Figure 5 <figure_flat_source_regions>`.

.. _figure_flat_source_regions:

.. figure:: ../../img/flat-source-regions.png
   :align: center
   :figclass: align-center
   :width: 400px

   **Figure 5**: A 4 :math:`\times` 4 lattice color-coded by flat source region.

.. note:: The runtime required by the plotting routine scales with the number of pixels in the image (the square of the ``gridsize`` parameter).
.. note:: The routine randomly selects a colormap at runtime. As a result, the colors in the figure will vary from run to run.

------------------
Flux Visualization
------------------

To plot the flat source region scalar fluxes, use the ``plotFluxes(...)`` routine in the ``openmoc.plotter`` module. The parameters accepted by this routine are described in :ref:`Table 10 <table_plot_fluxes>`.

.. _table_plot_fluxes:

=================  ===================  =========  =========  ============================================
Parameter          Type                 Default    Optional   Note
=================  ===================  =========  =========  ============================================
``geometry``       ``Geometry`` object  None       No         The geometry of interest
``solver``         ``Solver`` object    None       No         The solver used to converge the source
``energy_groups``  list                 [0]        No         Create separate plots for each energy group
``gridsize``       integer              250        Yes        The pixel resolution
=================  ===================  =========  =========  ============================================

**Table 10**: Parameters for the ``openmoc.plotter.plotFluxes(...)`` routine.

The code snippet below illustrates one possible configuration of parameters to the routine.

.. code-block:: python

    import openmoc.plotter as plot

    # Setup geomery and generate tracks
    ...

    # Setup solver and converge the source
    ...

    # Plot the fluxes for energy groups 1 and 7 in 500 x 500 pixel images
    plot.plotFluxes(geometry, solver, energy_groups=[1,7], gridsize=500)

A depiction of the group 1 and 7 fluxes for the C5G7 benchmark (:file:`/OpenMOC/sample-input/benchmarks/c5g7/c5g7.py`) is illustrated in :ref:`Figure 6 <figure_fluxes>`.

.. _figure_fluxes:

.. table:: 

   +------------------------------------------+-----------------------------------------+
   | .. _figa:                                | .. _figb:                               |
   |                                          |                                         |
   | .. image:: ../../img/flux-group-1.png    | .. image:: ../../img/flux-group-7.png   |
   |   :width: 72 %                           |   :width: 75 %                          |
   |   :align: center                         |   :align: center                        |
   +------------------------------------------+-----------------------------------------+ 


**Figure 6**: The fast and thermal fluxes in the C5G7 benchmark problem.


.. note:: The runtime required by the plotting routine scales with the number of pixels in the image (the square of the ``gridsize`` parameter).



.. _dictionary: http://docs.python.org/2/library/stdtypes.html#mapping-types-dict
.. _pickle: http://docs.python.org/2/library/pickle.html
.. _segmentation fault: http://en.wikipedia.org/wiki/Segmentation_fault
