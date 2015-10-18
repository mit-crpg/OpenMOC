##
# @file krylov.py
# @package openmoc.krylov
# @brief The krylov module provides the classes and routines for Krylov
#        subspace-based solution techniques.
# @details This module requires the SciPy Python package.
# @author William Boyd (wboyd@mit.edu)
# @date June 30, 2015

import sys

import numpy as np

import openmoc

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf


##
# @class krylov.py 'openmoc/krylov.py'
# @brief A Solver which uses a Krylov subspace-based method to solve for an
#        an arbitrary number of eigenmodes in a criticality problem.
# @details This class uses the Implicitly Restared Arnoldi Method (IRAM)
#        to compute the n highest order eigenvalues/vectors for a k-eigenvalue
#        criticality problem. This functionality is based on original work by
#        Colin Josey (cjosey@mit.edu).
#        NOTE: This functionality only works for vacuum boundary conditions.
#        FIXME: include all references used to construct this
class IRAMSolver(object):

  ##
  # @brief IRAMSolver class constructor
  # @param moc_solver an initialized OpenMOC Solver subclass (e.g. CPUSolver)
  def __init__(self, moc_solver):

    if not 'Solver' in type(moc_solver).__name__:
      py_printf('ERROR', 'Unable to initialize an IRAMSolver with %s ' + \
                'which is not an OpenMOC Solver object', str(moc_solver))
    
    self._moc_solver = moc_solver

    # Determine the floating point precision for Solver
    if self._moc_solver.isUsingDoublePrecision():
      self._precision = np.float64
    else:
      self._precision = np.float32

    # Determine if the user passed in a CUDA-enabled GPUSolver
    if 'GPUSolver' in type(moc_solver).__name__:
      self._with_cuda = True
    else:
      self._with_cuda = False

    # Compute the size of the LinearOperators used in the eigenvalue problem
    geometry = self._moc_solver.getGeometry()
    num_FSRs = geometry.getNumFSRs()
    num_groups = geometry.getNumEnergyGroups()
    self._op_size = num_FSRs * num_groups

    # Initialize solution-dependent class attributes to None
    self._num_modes = None
    self._interval = None
    self._outer_tol = None
    self._inner_tol = None
    self._A_op = None
    self._M_op = None
    self._F_op = None
    self._a_count = None
    self._m_count = None
    self._eigenvalues = None
    self._eigenvectors = None


  ##
  # @brief Compute all eigenmodes in the problem.
  # @details Uses the scipy.linalg package.
  # @param solver_mode the type of eigenmodes (FORWARD or ADJOINT) 
  # @param num_modes number of eigenmodes to compute
  # @param inner_method Krylov subspace method used for the Ax=b solve
  # @param outer_tol tolerance on the outer eigenvalue solve
  # @param inner_tol tolerance on the inner Ax=b solve
  # @param interval inner iteration interval for logging messages
  def computeEigenmodes(self, solver_mode=openmoc.FORWARD, num_modes=5, 
                        inner_method='gmres', outer_tol=1e-5, 
                        inner_tol=1e-6, interval=10):

    # Ensure that vacuum boundary conditions are used
    geometry = self._moc_solver.getGeometry()
    if (geometry.getMinXBoundaryType() != openmoc.VACUUM or
        geometry.getMaxXBoundaryType() != openmoc.VACUUM or
        geometry.getMinYBoundaryType() != openmoc.VACUUM or
        geometry.getMaxYBoundaryType() != openmoc.VACUUM):
        py_printf('ERROR', 'All boundary conditions must be ' + \
                           'VACUUM for the IRAMSolver')

    import scipy.sparse.linalg as linalg

    # Set solution-dependent class attributes based on parameters
    # These are accessed and used by the LinearOperators
    self._num_modes = num_modes
    self._inner_method = inner_method
    self._outer_tol = outer_tol
    self._inner_tol = inner_tol
    self._interval = interval

    # Initialize inner/outer iteration counters to zero
    self._m_count = 0
    self._a_count = 0

    # Initialize MOC solver
    self._moc_solver.initializePolarQuadrature()
    self._moc_solver.initializeExpEvaluator()
    self._moc_solver.initializeMaterials(solver_mode)
    self._moc_solver.initializeFluxArrays()
    self._moc_solver.initializeSourceArrays()
    self._moc_solver.initializeFSRs()
    self._moc_solver.countFissionableFSRs()
    self._moc_solver.zeroTrackFluxes()

    # Initialize SciPy operators
    op_shape = (self._op_size, self._op_size)
    self._A_op = linalg.LinearOperator(op_shape, self._A, dtype=self._precision)
    self._M_op = linalg.LinearOperator(op_shape, self._M, dtype=self._precision)
    self._F_op = linalg.LinearOperator(op_shape, self._F, dtype=self._precision)

    # Solve the eigenvalue problem
    timer = openmoc.Timer()
    timer.startTimer()
    vals, vecs = linalg.eigs(self._F_op, k=self._num_modes, tol=self._outer_tol)
    timer.stopTimer()

    # Print a timer report
    tot_time = timer.getTime()
    time_per_mode = tot_time / self._num_modes
    tot_time = '{0:.4e} sec'.format(tot_time)
    time_per_mode = '{0:.4e} sec'.format(time_per_mode)
    py_printf('RESULT', 'Total time to solution'.ljust(53, '.') + tot_time)
    py_printf('RESULT', 'Solution time per mode'.ljust(53, '.') + time_per_mode)

    # Store the eigenvalues and eigenvectors
    self._eigenvalues = vals
    self._eigenvectors = vecs

    # Restore the material data
    self._moc_solver.resetMaterials(solver_mode)



  ##
  # @brief Private routine for inner Ax=b solves with the scattering source.
  # @details Solves a fixed source problem using the scatter source for a given
  #          flux distribution. This corresponds to the left hand side of the
  #          generalized kAX = MX eigenvalue problem.
  # @param flux the flux used to compute the scattering source
  # @return the residual between input and computed fluxes
  def _A(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)
    flux_old = np.copy(flux)
    
    # Apply operator to flux
    self._a_count += 1
    self._moc_solver.setFluxes(flux)
    self._moc_solver.scatterTransportSweep()
    flux = self._moc_solver.getFluxes(self._op_size)

    # Print report to screen to update user on progress
    if self._a_count % self._interval == 0:
      py_printf('NORMAL', "Performed A operator sweep number %d", self._a_count)
    else:
      py_printf('INFO', "Performed A operator sweep number %d", self._a_count)
      
    # Return flux residual
    return flux_old - flux


  ##
  # @brief Private routine for inner Ax=b solves with the fission source.
  # @details Solves a fixed source problem using the fission source for a given
  #          flux distribution. This corresponds to the right hand side of the
  #          generalized kAX = MX eigenvalue problem.
  # @param flux the flux used to compute the fission source
  # @return the new flux computed from the MOC transport sweep
  def _M(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)

    # Apply operator to flux
    self._m_count += 1
    self._moc_solver.setFluxes(flux)
    self._moc_solver.fissionTransportSweep()
    flux = self._moc_solver.getFluxes(self._op_size)

    py_printf('NORMAL', "Performed M operator sweep number %d", self._m_count)

    # Return new flux
    return flux


  ##
  # @brief Private routine for outer eigenvalue solver method.
  # @details Uses a Krylov subspace method (e.g., GMRES, BICGSTAB) from the
  #          scipy.linalg package to solve the AX=B fixed scatter source
  #          problem.
  # @param flux a flux vector from the scipy.linalg.eigs routine
  # @return the flux computed from fission/scatter fixed source calculations
  def _F(self, flux):

    import scipy.sparse.linalg as linalg

    # Apply operator to flux - get updated flux from fission source
    flux = self._M_op * flux

    # Solve AX=B fixed scatter source problem using Krylov subspace method
    if self._inner_method == 'gmres':
      flux, x = linalg.gmres(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'lgmres':
      flux, x = linalg.lgmres(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'bicgstab':
      flux, x = linalg.bicgstab(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'cgs':
      flux, x = linalg.cgs(self._A_op, flux, tol=self._inner_tol)
    else:
      py_printf('ERROR', 'Unable to use %s to solve Ax=b', self._inner_method)

    # Check that solve completed without error before returning new flux
    if x != 0:
      py_printf('ERROR', 'Unable to solve Ax=b with %s', self._inner_method)
    else:
      return flux
