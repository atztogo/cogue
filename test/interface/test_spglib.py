import unittest

from cogue.interface.vasp_io import read_poscar_yaml
from cogue.interface.spglib import get_crystallographic_cell, get_primitive_cell

class TestSpglib(unittest.TestCase):

    def setUp(self):
        filename = "POSCAR.yaml"
        self._cell, poscar_order = read_poscar_yaml(filename)
    
    def tearDown(self):
        pass
    
    def test_get_crystallographic_cell(self):
        cell = get_crystallographic_cell(self._cell)
        print cell

    def test_get_primitive_cell(self):
        cell = get_primitive_cell(self._cell)
        print cell

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSpglib)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
