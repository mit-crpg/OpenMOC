import unittest
import numpy

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc

class TestMaterials(unittest.TestCase):

    def setUp(self):
        sigma_f = numpy.array([0.000625, 0.135416667])
        nu_sigma_f = numpy.array([0.0015, 0.325])
        sigma_s = numpy.array([[0.1, 0.117], [0., 1.42]])
        chi = numpy.array([1.0, 0.0])
        sigma_t = numpy.array([0.2208, 1.604])

        self.test_material = openmoc.Material()
        self.test_material.setName('2-group infinite medium')
        self.test_material.setNumEnergyGroups(2)
        self.test_material.setSigmaF(sigma_f)
        self.test_material.setNuSigmaF(nu_sigma_f)
        self.test_material.setSigmaS(sigma_s.flat)
        self.test_material.setChi(chi)
        self.test_material.setSigmaT(sigma_t)

    def test_fission_matrix(self):

          self.test_material.buildFissionMatrix()
          self.assertAlmostEqual(self.test_material.getFissionMatrixByGroup(1,1), .0015)

          self.assertAlmostEqual(self.test_material.getFissionMatrixByGroup(2,1), .325)

    def test_volume(self):

        self.test_material.setVolume(5)
        self.assertEqual(self.test_material.getVolume(), 5)
        
        self.test_material.incrementVolume(-2)
        self.assertEqual(self.test_material.getVolume(), 3)

        self.test_material.incrementVolume(2)
        self.assertEqual(self.test_material.getVolume(), 5)

if __name__ == '__main__':
    unittest.main()