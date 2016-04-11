#!/usr/bin/env python

import os
import sys
import math
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import LatticeGridInput

import openmoc

class QuadratureTestHarness(TestHarness):
    """Tests tracking over a lattice geometry."""

    def __init__(self):
        super(QuadratureTestHarness, self).__init__()
        self.quadratures = list()
        quadratures = openmoc.Quadrature.__subclasses__()
        for quadrature in quadratures:
          self.quadratures.append(quadrature())
        self._result = ''

    def _setup(self):
      return

    def _run_openmoc(self):
        """Segment tracks over the geometry and save the result to a string"""

        # Segmentize tracks over the geometry
        for azim in [4, 8, 16, 32, 64]:
            for polar in [4, 6]:
                for quad in self.quadratures:
                    sum_weights = 0
                    int_sine = 0
                    quad.setNumAzimAngles(azim)
                    quad.setNumPolarAngles(polar)
                    quad.initialize()
                    for a in range(int(azim/4)):
                        quad.setAzimSpacing(1.0, a)
                    quad.precomputeWeights(False)
                    for a in range(int(azim/2)):
                        for p in range(polar):
                            sum_weights += quad.getWeight(a, p) / \
                                math.sin(quad.getTheta(a, p))
                            int_sine += quad.getWeight(a, p)
                    if abs(sum_weights - 4 * math.pi) > 0.005:
                        self._result += 'Calculated ' + str(sum_weights) + \
                            ' for sum of weights for ' + type(quad).__name__\
                            + ' with ' + str(azim) + ' azimuthal angles and '\
                            + str(polar) + ' which exceeds tolerance of 0.005'\
                            ' from 4 PI\n'
                    if abs(int_sine - 9.8696) > 0.5:
                        self._result += 'Calculated ' + str(int_sine)\
                            + ' for the integral of sine for '\
                            + type(quad).__name__ + ' with ' + str(azim)\
                            + ' azimuthal angles and ' + str(polar)\
                            + ' which exceeds tolerance of 0.5 from 9.8696\n'

        if (self._result == ''):
            self._result += 'All Quadrature sets correctly summed their ' \
                'weights to 4 PI with tolerance 0.005 and calculated the ' \
                'integral of sine as 9.8696 with tolerance 0.5'



    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=True, num_segments=True, num_tracks=True,
                     hash_output=False):
        """Return the result string"""
        return self._result


if __name__ == '__main__':
    harness = QuadratureTestHarness()
    harness.main()
