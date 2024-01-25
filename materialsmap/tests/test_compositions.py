import unittest
from materialsmap.core.compositions import generateCompositions,createComposition

class TestCalculations(unittest.TestCase):

    def test_generateCompositions(self):
        composition_lists = generateCompositions(['Ag', 'Al'],1)
        self.assertEqual(composition_lists, [{'Ag': 0.0, 'Al': 0.0}, {'Ag': 0.0, 'Al': 1.0}, {'Ag': 1.0, 'Al': 0.0}],'The generateCompositions is wrong.')
    
    def test_createComposition(self):
        mater = {}
        mater['Cu'] = {'Cu':1}
        mater['Ag'] = {'Ag':1}
        mater['Al'] = {'Al':1}
        Compositions, numPoint, comp, numSimultion = createComposition(['Ag', 'Al'],['Cu', 'Ag', 'Al'],[{'Ag': 0.0, 'Al': 0.0}, {'Ag': 0.0, 'Al': 1.0}, {'Ag': 1.0, 'Al': 0.0}],mater,'')
        with self.subTest():
            self.assertEqual(Compositions,{'AL': [0, 1.0, 0], 'CU': [1.0, 0, 0], 'AG': [0, 0, 1.0], 'alloy_Cu': [1.0, 0.0, 0.0], 'alloy_Ag': [0.0, 0.0, 1.0], 'alloy_Al': [0.0, 1.0, 0.0], 'Index': [0, 1, 2]})
        with self.subTest():
            self.assertEqual(numPoint,3)
        with self.subTest():
            self.assertEqual(sorted(comp),['AG', 'AL', 'CU'])       
        with self.subTest():
            self.assertEqual(numSimultion,3)               


if __name__ == '__main__':
    unittest.main()