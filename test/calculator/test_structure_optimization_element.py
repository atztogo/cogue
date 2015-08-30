import unittest

import cogue.calculator.vasp as vasp
import numpy as np

class TestVasp(unittest.TestCase):

    def setUp(self):
        self._task = vasp.StructureOptimizationElement()
   
    def tearDown(self):
        pass
    
    def test_collect(self):
        for t in self._task:
            pass
        print self._task.get_current_cell()

if __name__ == '__main__':
    unittest.main()
