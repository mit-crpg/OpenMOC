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

   # Instantiate an OpenMOC Material class object
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
the geometry model, each unique closed volume in defined by its bounding
surfaces. The CSG formulation used in OpenMOC is described in more detail in :ref:`Constructive Solid Geometry <constructive_solid_geometry>`.

Cells
-----

Universes
---------

Lattices
--------



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
