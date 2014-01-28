.. _usersguide_input:

==========================
Writing Python Input Files
==========================

OpenMOC is provided to users as a Python API. As a result, there are not strict constraints on how an input file is written for an OpenMOC simulation as there are in many other scientific simulation codes. Instead, users may write a Python script or program and import OpenMOC and simply use the classes or routines which are necessary for a particular simulation. The :file:`OpenMOC/sample-input/` directory in the OpenMOC folder includes many example scripts for simulations ranging in complexity from a simple `pin cell`_ to the `C5G7 benchmark problem`_.

The following sections describe the essential portions of the OpenMOC API needed for reactor eigenvalue calculations. For more detail on the full extent of OpenMOC capabilities, users should reference the :ref:`OpenMOC API documentation <api>`.

.. note:: It is highly suggested that users acquire a basic understanding of Python before developing OpenMOC simulations. For users familiar with basic programming constructs such as loops and conditionals, the official `Python Tutorial`_ is an excellent place to learn Python basics. For users new to programming, the `Code Academy Python Course`_ provides an introduction to both programming essentials and the Python language.


-----------------------
Materials Specification
-----------------------

OpenMOC uses multi-group nuclear cross-sections prepared by some upstream processing tool such as the NJOY_ GROUPR module. In OpenMOC, cross-section data is encapsulated by the ``Material`` class in the main ``openmoc`` Python module. A ``Material`` class may be instantiated in Python and cross-sections may be loaded into it using NumPy_ data arrays as illustrated by the following code snippet:

.. code-block:: python

   import openmoc
   import numpy

   # Initialize material cross-sections using NumPy data arrays
   num_groups = 8
   sigma_a = numpy.array([0.1,0.15,0.2,0.25,0.35,0.4,0.45,0.5])
   sigma_f = numpy.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
   ...

   # Instantiate an OpenMOC Material class object with an
   # automatically-generated unique ID
   material = openmoc.Material(openmoc.material_id())

   # Set the number of energy groups in the material
   material.setNumEnergyGroups(num_groups)

   # Load the cross-section data into the material
   material.setSigmaA(sigma_a)
   material.setSigmaT(sigma_f)
   ...


For many simulations, defining the nuclear data cross-sections by hand in a Python script is cumbersome and error-prone. As a result, OpenMOC includes the ``openmoc.materialize`` module for importing nuclear data cross-sections from an HDF5_ or a Python pickle_ binary file. The ``materialize(...)`` routine is used to import data and instantiate ``Material`` objects returned via a Python dictionary_. The use of the ``openmoc.materialize`` module to import HDF5 and pickle binary files is illusrated in the following snippet:

.. code-block:: python

    import openmoc
    import openmoc.materialize as mat
    
    # Import cross-section data from an HDF5 file. This instantiates 
    # objects for each material and returns them in a dictionary
    # indexed by a name string defined in the pickle file.
    hdf5_materials = mat.materialize('materials-data.h5')

    # Retrieve the material called 'moderator' in the HDF5 file
    moderator = hdf5_materials['moderator']

    # Import cross-section data from a pickle file. This instantiates 
    # objects for each material and returns them in a dictionary
    # indexed by a name string defined in the pickle file
    pickle_materials = mat.materialize('materials-data.pkl')

    # Retrieve the material called 'fuel' in the pickle file
    fuel = pickle_materials['fuel']

The ``openmoc.materialize`` module defines a standard for cross-section data stored in binary files. First, each HDF5 file must end with the '.h5' or '.hdf5' extension. HDF5 files must include an `Energy Groups` attribute with the integer number of groups in the top level of the file data hierarchy. Finally, each material is defined as an `HDF5 group`_ with a string name to identify the material. Finally, the material group must contain the following floating point `HDF5 datasets`_:

  - 'Total XS'
  - 'Absorption XS'
  - 'Scattering XS'
  - 'Fission XS'
  - 'Nu Fission XS'
  - 'Chi'

The following code snippet illustrates the use of the h5py_ Python HDF5 interface to write an HDF5 file with material cross-section data adhering to the standard expected by the ``openmoc.materialize`` module:

.. code-block:: python

   import numpy
   import h5py

   # Create an HDF5 file to store multi-groups cross-sections
   f = h5py.File('materials-data.h5')

   # Set the number of energy groups
   f.attrs['Energy Groups'] = 8

   # Material 1

   # Create an HDF5 group for this material
   material_group = f.create_group('Material 1')

   # Initialize cross-sections as NumPy data arrays
   sigma_a = numpy.array([0.1,0.15,0.2,0.25,0.35,0.4,0.45,0.5])
   sigma_f = numpy.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
   ...

   # Create datasets for each cross-section type
   material_group.create_dataset('Absorption XS', data=sigma_a)
   material_group.create_dataset('Fission XS', data=sigma_f)
   ...

   # Material 2
   ...

   # Close and save the HDF5 file
   f.close()

Alternatively, for machine withouts HDF5 and/or h5py, materials data may be imported from a pickle_ binary file using the ``openmoc.materialize`` module. For pickle files, the materials data should be stored as a Python dictionary_. The dictionary must contain a key/value pair for the number of energy groups, and sub-dictionaries for each material's cross-sections. The following code snippet illustrates how one might populate a pickle file with material cross-section data adhering to the standard expected by the ``openmoc.materialize`` module:

.. code-block:: python

   import numpy
   import pickle

   # Initialize a Python dictionary to store the materials data
   data = dict()

   # Set the number of energy groups
   data['Energy Groups'] = 8

   # Material 1

   # Create a sub-dictoinary for this material
   data['Material 1'] = dict()

   # Initialize cross-sections as NumPy data arrays
   sigma_a = numpy.array([0.1,0.15,0.2,0.25,0.35,0.4,0.45,0.5])
   sigma_f = numpy.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
   ...

   # Create datasets for each cross-section type
   data['Material 1']['Absorption XS'] = sigma_a
   data['Material 1']['Fission XS'] = sigma_f
   ...

   # Material 2
   ...

   # Dump the Python dictionary of materials data to a pickle file
   pickle.dump(data, open('materials-data.pkl', 'wb'))


.. note:: Users must ensure that the total cross-section is equal to the absorption and scattering cross-section in each group. OpenMOC will throw a runtime error will be thrown if this condition does not hold true when materials are added to the ``Geometry`` object.


----------------------
Geometry Specification
----------------------

The geometry in OpenMOC is described using constructive solid geometry (CSG_),
also sometimes referred to as combinatorial geometry. CSG allows a user to
create complex objects using Boolean operators on a set of simpler surfaces. In
the geometry model, each unique closed volume is defined by its bounding
surfaces. The CSG formulation used in OpenMOC is described in more detail in :ref:`Constructive Solid Geometry <constructive_solid_geometry>`.

The following sections detail how to create surfaces, cells, universes and lattices to construct a simple 4 :math:`\times` 4 pin cell lattice. 


Surfaces
--------

In most cases, the first step towards building a reactor geometry is to create the surfaces defining boundaries between distinct regions. The CSG formulation for surfaces in OpenMOC is described in detail in :ref:`Surfaces and Halfspaces <surfaces-halfspaces>`. For LWRs, the most typical surfaces needed to model 2D rectangular lattices are the ``Circle``, ``XPlane``, and ``YPlane`` classes. The following code snippet illustrates how to create a circle to represent a fuel pin and reflective boundary planes to surround a 4 :math:`\times` 4 lattice.

.. code-block:: python

    # Initialize circular fuel pin surface
    circle = openmoc.Circle(x=0.0, y=0.0, radius=0.45)

    # Initialize the planar surfaces bounding the entire geometry
    left = openmoc.XPlane(x=-2.52)
    right = openmoc.XPlane(x=2.52)
    bottom = openmoc.YPlane(y=-2.52)
    top = openmoc.YPlane(y=2.52)

    # Set the boundary conditions for the bounding planes
    left.setBoundaryType(REFLECTIVE)
    right.setBoundaryType(REFLECTIVE)
    bottom.setBoundaryType(REFLECTIVE)
    top.setBoundaryType(REFLECTIVE)


Cells and Universes
-------------------

The next step to create a geometry is to instantiate cells which represent unique geometric shapes and use them to construct universes. The CSG formulations for cells and universes in OpenMOC are discussed in further detail in :ref:`Cells <cells>` and :ref:`Universes <universes>`, respectively. OpenMOC provides the ``CellBasic`` class for cells which are filled by a material. The following code snippet illustrates how to create cells filled by the fuel and moderator materials in the universe with ID = 1. Next, the script adds the appropriate halfspace of the circle surface created in the preceding section to each cell.

.. code-block:: python

    # Retrieve the IDs for the fuel and moderator materials
    uo2_id = materials['Fuel'].getId()
    water_id = materials['Water'].getId()

    # Initialize the cells for the fuel pin and moderator
    fuel = openmoc.CellBasic(universe=1, material=uo2_id)
    moderator = openmoc.CellBasic(universe=1, material=water_id)

    # Add the circle surface to each cell
    fuel.addSurface(halfspace=-1, surface=circle)
    moderator.addSurface(halfspace=+1, surface=circle)

In addition to cells filled with materials, OpenMOC provides the ``CellFill`` class for cells which may be filled with universes. As a result, a geometry may be constructed of a hierarchy of nested cells/universes. A hierarchichal geometry permits a simple treatment of repeating geometric structures on multiple length scales (e.g., rectangular arrays of fuel pins and fuel assemblies). 

OpenMOC does not place a limit on the hierarchical depth - or number of nested universe levels - that a user may define in constructing a geometry. The only limitation is that at the top of the hierarchy, a cell must be used to encapsulate the entire geometry in the universe with ID = 0. The following code snippet illustrates the creation of a ``CellFill`` which is filled by universe 10 - the lattice constructed in the next section - and which is part of universe 0. Finally, the appropriate halfspaces for the planes defined in the preceding section are added to the cell to enforce boundaries on the portion of universe 10 relevant to the geometry.

.. code-block:: python

    # Initialize a cell filled by the lattice universe. This cell 
    # resides within universe 0 which is designated for the top
    # level nested universe in the geometry.
    pin_cell_array = openmoc.CellFill(universe=0, universe_fill=10)

    # Add the bounding planar surfaces to each the cell containing
    # universe 0
    pin_cell_array.addSurface(halfspace=+1, left)
    pin_cell_array.addSurface(halfsapce=-1, right)
    pin_cell_array.addSurface(halfspace=+1, bottom)
    pin_cell_array.addSurface(halfspace=-1, top)


Lattices
--------

Once the cells for the geometry have been created, OpenMOC's ``Lattice`` class may be used to represent repeating patterns of the cells on a rectangular array. The CSG formulation for lattices is described further in :ref:`Lattices <lattices>`. In OpenMOC, the ``Lattice`` class is a subclass of the ``Universe`` class. The following code snippet illustrates the creation of a 4 :math:`\times` 4 lattice with each lattice cell filled by the universe with ID = 1. The total width and height of the lattice are defined as parameters when the lattice is initialized. The lattice dimensions are used to define the rectangular region of interest centered at the origin of each universe filling each lattice cell.

.. code-block:: python

    # Initialize the lattice for the geometry 
    lattice = openmoc.Lattice(id=10, width_x=5.04, width_y=5.04)

    # Assign each lattice cell a universe ID
    lattice.setLatticeCells([[1, 1, 1, 1],
                             [1, 1, 1, 1],
                             [1, 1, 1, 1],
                             [1, 1, 1, 1]])


Geometry
--------

The final step in creating a geometry is to instantiate OpenMOC's ``Geometry`` class. The ``Geometry`` class encapsulates all materials, surfaces, cells, universes and lattices. The following code snippet illustrates the creation of the geometry and the registration of each material, cell and lattice constructed in the preceding sections. The last line of the script is called once all primitives have been registered and is used to traverse the CSG hierarchy and index the flat source regions in the geometry.

.. code-block:: python

    # Initialize an empty geometry object
    geometry = openmoc.Geometry()

    # Add materials to the geometry first
    geometry.addMaterial(materials['Fuel'])
    geometry.addMaterial(materials['Water'])

    # Next, add all cells to the geometry
    geometry.addCell(fuel)
    geometry.addCell(moderator)
    geometry.addCell(pin_cell_array)

    # Next, add all lattices to the geometry
    geometry.addLattice(lattice)

    # Next, initialize the flat source regions in the geometry after
    # all materials, cells, and lattices have been added to it
    geometry.initializeFlatSourceRegions()


----------------
Track Generation
----------------

Once the geometry has been initialized for a simulation, the next step is to perform ray tracing for track generation. The track generation process and algorithms in OpenMOC are described in more detail in :ref:`Track Generation <track_generation>`. This step requires the instantiation of a ``TrackGenerator`` object and a function call to generate the tracks as illustrated in the following code snippet.

.. code-block:: python

    # Initialize the track generator after the geometry has been
    # constructed. Use 64 azimuthal angles and 0.05 cm track spacing.
    track_generator = openmoc.TrackGenerator(geometry, num_azim=64, \
                                             spacing=0.05)
    
    # Generate tracks using ray tracing across the geometry
    track_generator.generateTracks()


--------------------
MOC Source Iteration
--------------------

One of OpenMOC's ``Solver`` subclasses may be initialized given the ``Geometry`` and ``TrackGenerator`` objects created in the preceding sections. The most commonly used subclasses for OpenMOC simulations are itemized below:

  * ``ThreadPrivateSolver`` - multi-core CPUs, less memory efficient, excellent parallel scaling
  * ``CPUSolver`` - multi-core CPUs, memory efficient, poor parallel scaling
  * ``GPUSolver`` - GPUs, 30-50:math:`\times` faster than CPUs

The following code snippet illustrates the instantiation of the ``ThreadPrivateSolver`` for multi-core CPUs. The code assigns runtime parameters to the solver and calls the ``convergeSource(...)`` routine to execute the :ref:`MOC Source Iteration Algorithm <figure-overall-iterative-scheme>`.

.. code-block:: python

    # Initialize a solver for the simulation and set the number of
    # threads and source convergence threshold
    solver = openmoc.ThreadPrivateSolver(geometry, track_generator)
    solver.setNumThreads(4)
    solver.setSourceConvergenceThreshold(1E-5)

    # Converge the source with up to a maximum of 1000 source iterations
    solver.convergeSource(1000)

    # Print a report of the time to solution
    solver.printTimerReport()


.. _CSG: http://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _Python Tutorial: http://docs.python.org/2/tutorial/
.. _Code Academy Python Course: http://www.codecademy.com/tracks/python
.. _pin cell: https://github.com/mit-crpg/OpenMOC/tree/master/sample-input/pin-cell
.. _C5G7 benchmark problem: https://github.com/mit-crpg/OpenMOC/tree/master/sample-input/benchmarks/c5g7
.. _NumPy: http://www.numpy.org/
.. _NJOY: http://t2.lanl.gov/nis/njoy/title.html
.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _pickle: http://docs.python.org/2/library/pickle.html
.. _dictionary: http://docs.python.org/2/tutorial/datastructures.html#dictionaries
.. _h5py: http://www.h5py.org/
.. _HDF5 group: http://www.hdfgroup.org/HDF5/doc/UG/UG_frame09Groups.html
.. _HDF5 datasets: http://www.hdfgroup.org/HDF5/doc/UG/10_Datasets.html
