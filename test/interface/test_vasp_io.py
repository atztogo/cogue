import unittest

from cogue.interface.vasp_io import write_poscar, read_poscar, read_poscar_yaml, VaspCell

class TestVASPIO(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def test_read_poscar_yaml(self):
        filename = "POSCAR.yaml"
        cell, poscar_order = read_poscar_yaml(filename)
        for s, p in zip(cell.get_symbols(), cell.get_points().T):
            print s, p
        print write_poscar(cell)
        print poscar_order

if __name__ == '__main__':
    unittest.main()
