import unittest
import numpy as np
from materialsmap.core.GenerateEqScript import createEqScript
from materialsmap.core.ReadEqResult import getEqdata
from materialsmap.core.GenerateScheilScript import createScheilScript
from materialsmap.core.ReadScheilResult import getScheilSolidPhase
from importlib_resources import files
import os
from pathlib import Path
import json
import shutil 


class BaseTestCase(unittest.TestCase):
    
    def assertIsFile(self, path):
        if not Path(path).resolve().is_file():
            raise AssertionError("Thermo_Calc scripts didn't generate: %s" % str(path))

    def assertIsEmpty(self,path) -> None:
        msg = '{0} is empty.'.format(path)
        msg_exist = 'File does not exist: {0}'.format(path)
        if not Path(path).resolve().is_file():
            raise AssertionError(msg_exist)
        f = open(path,'r')
        try:
            obj = json.load(f)
        except:
            obj = None
        f.close()
        if obj is None or len(obj) == 0:
            raise AssertionError(msg)




class TestCalculations(BaseTestCase):
    """test thermo_calc functions

    Args:
        unittest (_type_): _description_
    """
    
    def setUp(self):
        """test for setting up the pycalphad eq and scheil simulations
        
        """
        self.path = str(files('materialsmap').joinpath('tests/testsCaseFiles'))

    def test_createEqScript(self):
        createEqScript(self.path)
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Equilibrium Simulation/0_CU_0.TCM')
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Equilibrium Simulation/1_AL_0.TCM')
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Equilibrium Simulation/2_AG_0.TCM')

    def test_createScheilScript(self):
        createScheilScript(self.path,2000)
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Scheil Simulation/0_CU-_0.TCM')
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Scheil Simulation/1_AL-_0.TCM')
        with self.subTest():
            self.assertIsFile(self.path+'/Thermo-calc/Scheil Simulation/2_AG-_0.TCM')

    def test_getEqdata(self):
        getEqdata(self.path)
        self.assertIsEmpty(self.path+'/Thermo-calc/Equilibrium Simulation/Result/data_mole.json')

    def test_getScheilSolidPhase(self):
        getScheilSolidPhase(self.path)
        self.assertIsEmpty(self.path+'/Thermo-calc/Scheil Simulation/Result/data_mole.json')

    @classmethod
    def tearDownClass(self):
        self.path = str(files('materialsmap').joinpath('tests/testsCaseFiles'))
        try:
            os.remove(self.path+'/Thermo-calc/Equilibrium Simulation/0_CU_0.TCM')
            os.remove(self.path+'/Thermo-calc/Equilibrium Simulation/1_AL_0.TCM')
            os.remove(self.path+'/Thermo-calc/Equilibrium Simulation/2_AG_0.TCM')
            os.remove(self.path+'/Thermo-calc/Scheil Simulation/0_CU-_0.TCM')
            os.remove(self.path+'/Thermo-calc/Scheil Simulation/1_AL-_0.TCM')
            os.remove(self.path+'/Thermo-calc/Scheil Simulation/2_AG-_0.TCM')
            os.remove(self.path+'/Thermo-calc/Equilibrium Simulation/Result/data_mole.json')
            os.remove(self.path+'/Thermo-calc/Scheil Simulation/Result/data_mole.json')
            os.remove(self.path+'/Thermo-calc/Scheil Simulation/Result/ScheilResults.xlsx')
        except OSError as why:
            print(why)


if __name__ == '__main__':
    unittest.main()