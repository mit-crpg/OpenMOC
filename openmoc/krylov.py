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
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs, gmres


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
    self._a_count = None
    self._m_count = None
    self._interval = None
    self._tolerance = None

    self._A_operator = None
    self._M_operator = None
    self._F_operator = None


  def computeEigenmodes(self, num_modes=5, tol=1e-5, interval=10,
                        max_outer=None, max_inner=None):

    self._tolerance = tol
    self._interval = interval

    self._m_count = 0
    self._a_count = 0

    self._solver.initializeMemory()
    size = self._solver.getOperatorSize()

    self._A_operator = LinearOperator((size, size), matvec=self._A, dtype='float64')
    self._M_operator = LinearOperator((size, size), matvec=self._M, dtype='float64')
    self._F_operator = LinearOperator((size, size), matvec=self._F, dtype='float64')

    vals, vecs = eigs(self._F_operator, k=num_modes, tol=self._tolerance)

    print vals


  def _A(self, flux):

    # Do not pass imaginary numbers to SWIG
    flux = np.real(flux).astype(np.float64)
    flux_old = copy.copy(flux)
    
    self._solver.scatterTransportSweep(flux)
    
    self._a_count += 1

    if self._a_count % self._interval == 0:
      py_printf('NORMAL', "Performed A operator sweep number %d", self._a_count)
    else:
      py_printf('INFO', "Performed A operator sweep number %d", self._a_count)
    
    return flux_old - flux


  def _M(self, flux):

    # Do not pass imaginary numbers to SWIG
    flux = np.real(flux).astype(np.float64)
    
    self._solver.fissionTransportSweep(flux)
    
    self._m_count += 1

    py_printf('NORMAL', "Performed M operator sweep number %d", self._m_count)
    
    return flux


  def _F(self, flux):
  
    flux = self._M_operator * flux

    # Note, gmres must converge beyond tolerance to work.
    flux, x = gmres(self._A_operator, flux, tol=self._tolerance/10) 
    
    return flux
