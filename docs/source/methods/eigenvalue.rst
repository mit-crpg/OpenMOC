.. _methods_eigenvalue:

=======================
Eigenvalue Calculations
=======================

An eigenvalue calculation, also referred to as a criticality calculation, is a
transport simulation wherein the source of neutrons includes a fissionable
material. Some common eigenvalue calculations include the simulation of nuclear
reactors, spent fuel pools, nuclear weapons, and other fissile systems. The
reason they are called *eigenvalue* calculations is that the transport equation
becomes an eigenvalue equation if a fissionable source is present since then the
source of neutrons will depend on the flux of neutrons itself.

This section will explore the theory behind and implementation of eigenvalue
calculations in a method of characteristics code.
