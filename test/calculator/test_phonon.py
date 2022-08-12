import unittest

import numpy as np

import cogue
import cogue.calculator.vasp as vasp

symbols = ["Si"] * 2 + ["O"] * 4
lattice = [[4.65, 0, 0], [0, 4.75, 0], [0, 0, 3.25]]
points = np.transpose(
    [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
)
cell = cogue.cell(lattice=lattice, points=points, symbols=symbols)


class TestVasp(unittest.TestCase):
    def setUp(self):
        self._task = vasp.electronic_structure(cell=cell)

    def tearDown(self):
        pass

    def test_collect(self):
        for t in self._task:
            pass


if __name__ == "__main__":
    unittest.main()
