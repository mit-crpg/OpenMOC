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

    def test_ids(self):

        openmoc.reset_material_id()
        material_2 = openmoc.Cell()
        self.assertEqual(material_2.getId(), 1000000)
        openmoc.maximize_material_id(10000000)
        material_3 = openmoc.Cell()
        self.assertEqual(material_3.getId(), 1000001)
        self.test_material.printString()

    def test_instances(self):

        self.assertEqual(self.test_material.getNumInstances(), 0)
        self.test_material.setNumInstances(99)
        self.assertEqual(self.test_material.getNumInstances(), 99)
        self.test_material.incrementNumInstances()
        self.assertEqual(self.test_material.getNumInstances(), 100)

    def test_fission_matrix(self):

        material = openmoc.Material()
        with self.assertRaises(Exception): material.buildFissionMatrix()

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

    def test_get_cross_section_by_group(self):

        self.assertAlmostEqual(self.test_material.getSigmaTByGroup(1), 0.2208)
        self.assertAlmostEqual(self.test_material.getSigmaAByGroup(1), 0.0038)
        self.assertAlmostEqual(self.test_material.getSigmaFByGroup(1), 0.000625)
        self.assertAlmostEqual(self.test_material.getSigmaSByGroup(1, 1), 0.1)
        self.assertAlmostEqual(self.test_material.getChiByGroup(1), 1.)
        self.assertAlmostEqual(self.test_material.getNuSigmaFByGroup(1), 0.0015)

    def test_get_cross_section(self):

        material = openmoc.Material()
        with self.assertRaises(Exception): material.getSigmaS()
        with self.assertRaises(Exception): material.getSigmaA()
        with self.assertRaises(Exception): material.getSigmaT()
        with self.assertRaises(Exception): material.getSigmaF()
        with self.assertRaises(Exception): material.getNuSigmaF()

    def test_cross_section_alignment(self):

        material = openmoc.Material()
        with self.assertRaises(Exception): material.alignData()

        self.test_material.alignData()
        self.assertEqual(self.test_material.isDataAligned(), True)
        self.assertEqual(self.test_material.getNumVectorGroups(), 1)

        self.assertAlmostEqual(self.test_material.getSigmaTByGroup(1), 0.2208)
        #FIXME SigmaA is not copied during alignment process
        self.assertAlmostEqual(self.test_material.getSigmaAByGroup(1), 0.1208)
        self.assertAlmostEqual(self.test_material.getSigmaFByGroup(1), 0.000625)
        self.assertAlmostEqual(self.test_material.getSigmaSByGroup(1, 1), 0.1)
        self.assertAlmostEqual(self.test_material.getChiByGroup(1), 1.)
        self.assertAlmostEqual(self.test_material.getNuSigmaFByGroup(1), 0.0015)

if __name__ == '__main__':
    unittest.main()
