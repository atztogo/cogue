import unittest

import cogue.calculator.vasp as vasp


class TestVasp(unittest.TestCase):
    def setUp(self):
        self._task = vasp.ElectronicStructure()
        # self._task = vasp.electronic_structure(cell=cell)

    def tearDown(self):
        pass

    def test_collect(self):
        for t in self._task:
            pass
        print(self._task.get_properties()["forces"])


if __name__ == "__main__":
    unittest.main()
