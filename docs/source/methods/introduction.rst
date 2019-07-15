.. _methods_introduction:

============
Introduction
============

The method of characteristics is one technique for solving partial differential equations.  MOC is one of the most common methods with real world applications in production lattice physics tools used today [Smith]_. The prospects for the MOC algorithm as an eventual successor to low-order diffusion-based method for reactor analysis are promising, and extensive research efforts into the efficient use of the algorithm are ongoing [Wu]_, [Taylor]_, [Liu]_, [Talamo]_. As a result, there is an opportunity for advanced parallel algorithms for high performance computing machines and nonlinear acceleration schemes to propel the application of MOC to full-core reactor physics calculations.

:ref:`Section 2 <method_of_characteristics>` describes how the method of characteristics is applied to solve the steady-state neutron transport equation in OpenMOC. :ref:`Section 3 <eigenvalue_calculations>` presents the algorithms implemented in OpenMOC to solve the MOC equations. :ref:`Section 4 <constructive_solid_geometry>` discusses the constructive solid geometry formulation used to represent geometric models in OpenMOC. :ref:`Section 5 <track_generation>` reviews the angular quadrature and track generation algorithm in OpenMOC. :ref:`Section 6 <parallelization>` outlines the parallel algorithms in OpenMOC for high performance computing platforms while :ref:`Section 7 <cmfd>` describes CMFD, a nonlinear acceleration method in OpenMOC. :ref:`Section 7 <krylov>` gives elements for replacing the fission iteration algorithm by Krylov subspace algorithms.


References
==========

.. [Smith] K. Smith and J. D. Rhodes, "CASMO-4 Characteristics Methods for Two-Dimensional PWR and BWR Core Calculations." *Transactions of the American Nuclear Society*, **83**, pp. 294 (2000).

.. [Wu] G. Wu and R. Roy, "A New Characteristics Algorithm for 3D Transport Calculations." *Annals of Nuclear Energy*, **30**, pp. 1-6 (2003).

.. [Taylor] J. Taylor, D. Knott and A. J. Baratta, "A Method of Characteristics Solution to the OECD/NEA 3D Neutron Transport Benchmark Problem." *Proceedings of the Joint International Topical Meeting on Mathematics and Computation and Supercomputing in Nuclear Applications*, Monterey, CA, USA (2007).

.. [Liu] Z. Liu, H. Wu, L. Cao, Q. Chen, Y. Li, "A New Three-Dimensional Method of Characteristics for the Neutron Tansport Calculation." *Annals of Nuclear Energy*, **38**, pp. 447-454 (2011).

.. [Talamo] A. Talamo, "Numerical Solution of the Time Dependent Neutron Transport Equation by the Method of the Characteristics." *Journal of Computational Physics}*, **240**, pp. 248-267 (2013).
