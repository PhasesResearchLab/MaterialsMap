import unittest
import os
import numpy as np
from fmap.core.GenerateEqScript import createEqScript
from fmap.core.ReadEqResult import getEqdata
from fmap.core.GenerateScheilScript import createScheilScript
from fmap.core.ReadScheilResult import getScheilSolidPhase

class TestCalculations(unittest.TestCase):
    
    def setUp(self):
        #path = resources.files('fmap').joinpath('tests/testsCaseFiles')
        path = './testsCaseFiles'
        settings= ((600, 2000, 10),3,3,['AL', 'CU', 'AG'],['Cu', 'Ag', 'Al'],['Ag', 'Al'],path+'/Ag-Al-Cu.TDB',101325,'massFraction')
        np.save(f'{path}/setting.npy',settings)
    def test_createEqScript(self):
        path = './testsCaseFiles'
        createEqScript(path)
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Equilibrium Simulation/0_AL_0.TCM')
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Equilibrium Simulation/1_AG_0.TCM')
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Equilibrium Simulation/2_CU_0.TCM')
    def test_createScheilScript(self):
        path = './testsCaseFiles'
        createScheilScript(path,2000)
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Scheil Simulation/0_AL-_0.TCM')
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Scheil Simulation/1_AG-_0.TCM')
        with self.subTest():
            self.assertTrue(path+'/Thermo-calc/Scheil Simulation/2_CU-_0.TCM')
    def test_getEqdata(self):
        path = './testsCaseFiles'
        getEqdata(path)
        self.assertTrue(path+'/Thermo-calc/Equilibrium Simulation/Result/data.json')
    def test_getScheilSolidPhase(self):
        path = './testsCaseFiles'
        getScheilSolidPhase(path)
        self.assertTrue(path+'/Thermo-calc/Scheil Simulation/Result/data.json')
    # def tearDown(self):
    #     path = './testsCaseFiles'
    #     os.remove(path+'/Thermo-calc')


if __name__ == '__main__':
    unittest.main()