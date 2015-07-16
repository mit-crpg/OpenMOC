##
# @file krylov.py
# @package openmoc.krylov
# @brief The krylov module provides the classes and routines for Krylov
#        subspace-based solution techniques.
# @details NOTE: This module requires the SciPy Python package.
# @author William Boyd (wboyd@mit.edu)
# @date June 30, 2015

import sys
import copy
import numpy as np
import scipy.sparse.linalg as linalg

import openmoc

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from log import *
# For Python 3.X.X
else:
  from openmoc.log import *


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
  # @param solver an initialized OpenMOC Solver subclass (e.g. CPUSolver)
  def __init__(self, solver):

    if not 'Solver' in type(solver).__name__:
      py_printf('ERROR', 'Unable to initialize an IRAMSolver with %s ' + \
                'which is not an OpenMOC Solver object', str(solver))
    
    self._solver = solver

    if self._solver.isUsingDoublePrecision():
      self._precision = np.float64
    else:
      self._precision = np.float32

    if 'GPUSolver' in type(solver).__name__:
      self._with_cuda = True
    else:
      self._with_cuda = False

    geometry = self._solver.getGeometry()
    num_FSRs = geometry.getNumFSRs()
    num_groups = geometry.getNumEnergyGroups()
    self._op_size = num_FSRs * num_groups

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
  # @param num_modes number of eigenmodes to compute
  # @param inner_method Krylov subspace method used for the Ax=b solve
  # @param outer_tol tolerance on the outer eigenvalue solve
  # @param inner_tol tolerance on the inner Ax=b solve
  # @param interval inner iteration interval for logging messages
  def computeEigenmodes(self, num_modes=5, inner_method='gmres', 
                        outer_tol=1e-5, inner_tol=1e-6, interval=10):

    self._inner_method = inner_method
    self._outer_tol = outer_tol
    self._inner_tol = inner_tol
    self._interval = interval

    self._m_count = 0
    self._a_count = 0

    # Initialize MOC solver
    self._solver.initializePolarQuadrature()
    self._solver.initializeExpEvaluator()
    self._solver.initializeFluxArrays()
    self._solver.initializeSourceArrays()
    self._solver.initializeFSRs()
    self._solver.countFissionableFSRs()
    self._solver.zeroTrackFluxes()

    # Initialize SciPy operators
    op_shape = (self._op_size, self._op_size)
    self._A_op = linalg.LinearOperator(op_shape, self._A, dtype=self._precision)
    self._M_op = linalg.LinearOperator(op_shape, self._M, dtype=self._precision)
    self._F_op = linalg.LinearOperator(op_shape, self._F, dtype=self._precision)

    # Solve the eigenvalue problem
    timer = openmoc.Timer()
    timer.startTimer()
    vals, vecs = linalg.eigs(self._F_op, k=num_modes, tol=self._outer_tol)
    timer.stopTimer()

    # Print a timer report
    tot_time = timer.getTime()
    time_per_mode = tot_time / num_modes
    tot_time = '{0:.4e}'.format(tot_time)
    time_per_mode = '{0:.4e}'.format(time_per_mode)
    py_printf('RESULT', 'Total time to solution'.ljust(53, '.')+tot_time)
    py_printf('RESULT', 'Solution time per mode'.ljust(53, '.')+time_per_mode)

    # Store the eigenvalues and eigenvectors
    self._eigenvalues = vals
    self._eigenvectors = vecs


  ##
  # @brief Private routine for inner Ax=b solves with the scattering source.
  # @details Applies a transport sweep to the scattering source for a given
  #          flux distribution. This corresponds to the left hand side of the
  #          generalized kAX = MX eigenvalue problem.
  # @param flux the flux used to compute the scattering source
  # @return the residual between input flux and output flux
  def _A(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)
    flux_old = np.copy(flux)
    
    # Apply operator to flux
    self._a_count += 1
    self._solver.setFluxes(flux)
    self._solver.scatterTransportSweep()
    flux = self._solver.getFluxes(self._op_size)

    # Print report to screen to update user on progress
    if self._a_count % self._interval == 0:
      py_printf('NORMAL', "Performed A operator sweep number %d", self._a_count)
    else:
      py_printf('INFO', "Performed A operator sweep number %d", self._a_count)
      
    # Return flux residual
    return flux_old - flux


  ##
  # @brief Private routine for inner Ax=b solves with the fission source.
  # @details Applies a transport sweep to the fission source for a given
  #          flux distribution. This corresponds to the right hand side of the
  #          generalized kAX = MX eigenvalue problem.
  # @param flux the flux used to compute the fission source
  # @return the new flux computed from the MOC transport sweep
  def _M(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)

    # Apply operator to flux
    self._m_count += 1
    self._solver.setFluxes(flux)
    self._solver.fissionTransportSweep()
    flux = self._solver.getFluxes(self._op_size)

    py_printf('NORMAL', "Performed M operator sweep number %d", self._m_count)

    # Return new flux
    return flux


  ##
  # @brief Private routine for outer eigenvalue 
  # @details Uses a Krylov subspace method (e.g., GMRES, BICGSTAB) to solve 
  #          the Ax=b problem.
  # @param flux 
  # @return the new flux
  def _F(self, flux):
  
    # Apply operator to flux
    flux = self._M_op * flux

    # Note, gmres must converge beyond tolerance to work.
    if self._inner_method == 'gmres':
      flux, x = linalg.gmres(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'lgmres':
      flux, x = linalg.lgmres(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'cg':
      flux, x = linalg.cg(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'bicgstab':
      flux, x = linalg.bicgstab(self._A_op, flux, tol=self._inner_tol)
    elif self._inner_method == 'cgs':
      flux, x = linalg.cgs(self._A_op, flux, tol=self._inner_tol)
    else:
      py_printf('ERROR', 'Unable to use %s to solve Ax=b', self._inner_method)

    # Check that Ax=b solve completed without error before returning new flux
    if x != 0:
      py_printf('ERROR', 'Unable to solve Ax=b with %s', self._inner_method)
    else:
      return flux
