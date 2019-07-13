.. _software_design:

===============
Software Design
===============

OpenMOC uses a compiled language coupled with a scripting language "glue" [Sanner]_, a methodology that has increasingly gained traction across disciplines since it enables both usability and performance. The majority of the source code is written in C/C++ as it is the most robust and well supported general purpose, high performance, compiled programming language with object-oriented features. In addition, OpenMOC's solver routines for the GPU are written in NVIDIA's CUDA programming language [CUDA]_ - a compiled language with similar syntax to C/C++. The widely adopted Simplified Wrapper Interface Generator (SWIG) [Beazley]_ is deployed to expose the C/C++/CUDA classes and routines to the Python scripting language. The model to couple the compiled C/C++/CUDA code to the Python scripting language is illustrated in :ref:`Figure 1 <figure_software_design>`. OpenMOC's Python interface allows for rapid prototyping and testing of code features and tight integration with the rich ecosystem of powerful data processing and visualization tools developed for Python.

.. _figure_software_design:

.. figure:: ../../img/openmoc-software-design.png
   :align: center
   :figclass: align-center
   :width: 500px

   **Figure 1**: Programming model in OpenMOC.

----------------------
Object-Oriented Design
----------------------

OpenMOC is designed using the Object-Oriented (OO) programming paradigm, the standard for software development for over two decades. In OO programming, data structures called **classes** are used to encapsulate both attributes and subroutines. Object-oriented programming is powerful since it directly enables code generalization and requires the use of **trust boundaries** which lead to more resilient code. An OpenMOC simulation is created through the instantiation and manipulation of class objects. A complete listing of classes in OpenMOC is displayed in :ref:`Table 1 <table_classes>`.

.. _table_classes:


===========================  =======================  ============================
Class                        Parent Class             Category
===========================  =======================  ============================
``Surface``                  N/A                      Constructive Solid Geometry
``Plane``                    ``Surface``              Constructive Solid Geometry
``XPlane``                   ``Plane``                Constructive Solid Geometry
``YPlane``                   ``Plane``                Constructive Solid Geometry
``ZCylinder``                ``Surface``              Constructive Solid Geometry
``Cell``                     N/A                      Constructive Solid Geometry
``Universe``                 ``Universe``             Constructive Solid Geometry
``Lattice``                  N/A                      Constructive Solid Geometry
``Geometry``                 N/A                      Constructive Solid Geometry
``TrackGenerator``           N/A                      Ray Tracing
``TrackGenerator3D``         TrackGenerator           Ray Tracing
``Quadrature``               N/A                      Ray Tracing
``Material``                 N/A                      Nuclear Data
``Solver``                   N/A                      Method of Characteristics
``CPUSolver``                ``Solver``               Method of Characteristics
``CPULSSolver``              ``CPUSolver``            Method of Characteristics
``GPUSolver``                ``Solver``               Method of Characteristics
``VectorizedSolver``         ``CPUSolver``            Method of Characteristics
===========================  =======================  ============================

**Table 1**: OpenMOC's classes.


--------------
Python Modules
--------------

OpenMOC’s Python interface makes it relatively easy to create complicated simulation input. Generating input for an OpenMOC simulation does not involve writing an input file in the traditional sense. OpenMOC leverages the flexiblity provided by Python to allow users complete control to build their inputs in one or more scripts just as one may do for any Python program. The user imports the necessary OpenMOC modules (see :ref:`Table 2 <table_python_modules>`) into Python and “builds” a simulation using only those classes and routines which are needed.


.. _table_python_modules:

=======================  ======================================
Module                   Description
=======================  ======================================
``openmoc``              The main module for OpenMOC
``openmoc.options``      Command line options
``openmoc.log``          Level-based logging messages
``openmoc.materialize``  Import/export multi-group nuclear data
``openmoc.plotter``      Visualizations for geometry, flux, etc.
``openmoc.process``      Data processing
``openmoc.cuda``         Solver for NVIDIA GPUs
``openmoc.krylov``       Solver based on krylov IRAM methods
=======================  ======================================

**Table 2**: OpenMOC's Python modules.


.. [Sanner] M. Sanner, "Python: A Programming Language for Software Integration and Development." Journal of Molecular Graphics and Modelling, **17(1)**, pp. 57-61 (1999).
.. [CUDA] NVIDIA, "NVIDIA CUDA C Programming Guide." http://docs.nvidia.com/cuda/cuda-c-programming-guide/ (2013).
.. [Beazley] D. Beazley, "Automated Scientific Software Scripting with SWIG." Future Generation Computer Systems, **19(5)**, pp. 599-609 (2003).
