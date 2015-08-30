import unittest

import numpy as np
from cogue.crystal.cell import Cell, sort_cell_by_symbols

class TestCell(unittest.TestCase):

    def setUp(self):
        symbols = ['Si', 'O', 'O', 'Si', 'O', 'O']
        lattice = [[4.65, 0, 0],
                   [0, 4.75, 0],
                   [0, 0, 3.25]]
        points=np.transpose([[0.0, 0.0, 0.0],
                             [0.3, 0.3, 0.0],
                             [0.7, 0.7, 0.0],
                             [0.5, 0.5, 0.5],
                             [0.2, 0.8, 0.5],
                             [0.8, 0.2, 0.5]])
        self._cell = Cell(lattice=lattice,
                          points=points,
                          symbols=symbols)
    
    def tearDown(self):
        pass
    
    def test_str(self):
        print self._cell

    def test_sort_cell_by_symbols(self):
        print sort_cell_by_symbols(self._cell)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCell)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
