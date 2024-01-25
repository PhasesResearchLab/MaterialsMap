import unittest
from pathlib import Path
import numpy as np
from materialsmap.core.pycalphad_run import pycalphad_eq,pycalphad_scheil
from importlib_resources import files
import json
import shutil 

class BaseTestCase(unittest.TestCase):
    
    def assertIsEmpty(self,path) -> None:
        msg = '{0} is not empty.'.format(path)
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
    """test pycalphad functions

    Args:
        unittest (_type_): _description_
    """
    
    def setUp(self):
        """test for setting up the pycalphad eq and scheil simulations
        
        """
        self.path = str(files('materialsmap').joinpath('tests/testsCaseFiles'))
        settings= ((600, 2000, 10),3,3,['AL', 'CU', 'AG'],['Cu', 'Ag', 'Al'],['Ag', 'Al'],self.path+'/Ag-Al-Cu.TDB',101325,'massFraction')
        np.save(self.path+'/setting.npy', settings)

    def test_pycalphad_eq(self):
        pycalphad_eq(self.path)
        self.assertIsEmpty(self.path+'/Pycalphad/Equilibrium Simulation/Result/data_mole.json')

    def test_pycalphad_scheil(self):
        pycalphad_scheil(self.path,2000) 
        self.assertIsEmpty(self.path+'/Pycalphad/Scheil Simulation/Result/data_mole.json')    

    @classmethod
    def tearDownClass(self):
        self.path = str(files('materialsmap').joinpath('tests/testsCaseFiles'))
        try:
            shutil.rmtree(self.path+'/Pycalphad')
        except OSError as why:
            print(why)

if __name__ == '__main__':
    unittest.main()