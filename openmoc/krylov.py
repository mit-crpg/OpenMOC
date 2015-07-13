##
# @file krylov.py
# @package openmoc.krylov
# @brief The krylov module provides the classes and routines for Krylov
#        subspace-based solution techniques.
# @details NOTE: This module requires the SciPy Python package.
# @author William Boyd (wboyd@mit.edu)
# @date June 30, 2015

import sys

## @var openmoc
#  @brief The openmoc module in use in the Python script using the
#         openmoc.krylov module.
openmoc = ''

# Determine which OpenMOC module is being used
if 'openmoc.gnu.double' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.double']
elif 'openmoc.gnu.single' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.single']
elif 'openmoc.intel.double' in sys.modules:
  openmoc = sys.modules['openmoc.intel.double']
elif 'openmoc.intel.single' in sys.modules:
  openmoc = sys.modules['openmoc.intel.single']
elif 'openmoc.bgq.double' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.double']
elif 'openmoc.bgq.single' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.single']
else:
  import openmoc

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from log import *
# For Python 3.X.X
else:
  from openmoc.log import *

import copy
import numpy as np
import scipy.sparse.linalg as linalg


# TODO: Use PyCuda to create handle on array shared between CPU/GPU?
# TODO: Remove CPUSolver::putFluxes(...) in place of storing array pointer
# TODO: Add a Timer report

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
  def __init__(self, solver):

    if not isinstance(solver, openmoc.Solver):
      py_printf('ERROR', 'Unable to initialize an IRAMSolver with %s ' + \
                'which is not an OpenMOC Solver object', str(solver))
    
    self._solver = solver

    if self._solver.isUsingDoublePrecision():
      self._precision = np.float64
    else:
      self._precision = np.float32

    self._size = self.getOperatorSize()

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
    op_shape = (self._size, self._size)
    self._A_op = linalg.LinearOperator(op_shape, self._A, dtype=self._precision)
    self._M_op = linalg.LinearOperator(op_shape, self._M, dtype=self._precision)
    self._F_op = linalg.LinearOperator(op_shape, self._F, dtype=self._precision)

    # Solve the eigenvalue problem
    openmoc.Timer.startTimer()
    vals, vecs = linalg.eigs(self._F_op, k=num_modes, tol=self._outer_tol)
    openmoc.Timer.stopTimer()
    openmoc.Timer.recordSplit('Total time')
    tot_time = openmoc.Timer.getTime('Total time')
    py_printf('RESULT', 'Total time to solution'.ljust(53, '.') + str(tot_time))

    self._eigenvalues = vals
    self._eigenvectors = vecs


  def getOperatorSize(self):
    geometry = self._solver.getGeometry()
    num_FSRs = geometry.getNumFSRs()
    num_groups = geometry.getNumEnergyGroups()
    return num_FSRs * num_groups


  def _A(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)
    flux_old = np.copy(flux)
    
    # Apply operator to flux
    self._a_count += 1
    self._solver.setFluxes(flux)
    self._solver.scatterTransportSweep()

    if self._a_count % self._interval == 0:
      py_printf('NORMAL', "Performed A operator sweep number %d", self._a_count)
    else:
      py_printf('INFO', "Performed A operator sweep number %d", self._a_count)
    
    return flux_old - flux


  def _M(self, flux):

    # Remove imaginary components from NumPy array
    flux = np.real(flux).astype(self._precision)
    
    # Apply operator to flux
    self._m_count += 1
    self._solver.setFluxes(flux)
    self._solver.fissionTransportSweep()

    py_printf('NORMAL', "Performed M operator sweep number %d", self._m_count)
    
    return flux


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

    if x != 0:
      py_printf('ERROR', 'Unable to solve Ax=b with %s', self._inner_method)
    else:
      return flux
