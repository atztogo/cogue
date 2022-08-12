import unittest

import numpy as np

from cogue.crystal.cell import Cell, sort_cell_by_symbols


class TestCell(unittest.TestCase):
    def setUp(self):
        symbols = ["Si", "O", "O", "Si", "O", "O"]
        lattice = [[4.65, 0, 0], [0, 4.75, 0], [0, 0, 3.25]]
        points = np.transpose(
            [
                [0.0, 0.0, 0.0],
                [0.3, 0.3, 0.0],
                [0.7, 0.7, 0.0],
                [0.5, 0.5, 0.5],
                [0.2, 0.8, 0.5],
                [0.8, 0.2, 0.5],
            ]
        )
        self._cell = Cell(lattice=lattice, points=points, symbols=symbols)

    def tearDown(self):
        pass

    def _test_str(self):
        print(self._cell)

    def _test_sort_cell_by_symbols(self):
        print(sort_cell_by_symbols(self._cell))

    def test_lattice(self):
        lattice = self._cell.lattice
        self._cell.lattice = lattice
        np.testing.assert_allclose(lattice, self._cell.lattice)
        self._cell.set_lattice(lattice)
        self._cell.get_lattice()
        np.testing.assert_allclose(lattice, self._cell.get_lattice())

    def test_points(self):
        points = self._cell.points
        self._cell.points = points
        np.testing.assert_allclose(points, self._cell.points)
        self._cell.set_points(points)
        self._cell.get_points()
        np.testing.assert_allclose(points, self._cell.get_points())

    def test_numbers(self):
        numbers = self._cell.numbers
        self._cell.numbers = numbers
        np.testing.assert_array_equal(numbers, self._cell.numbers)
        self._cell.set_numbers(numbers)
        self._cell.get_numbers()
        np.testing.assert_array_equal(numbers, self._cell.get_numbers())

    def test_symbols(self):
        symbols = self._cell.symbols
        self._cell.symbols = symbols
        for s, ss in zip(symbols, self._cell.symbols):
            self.assertEqual(s, ss)
        self._cell.set_symbols(symbols)
        self._cell.get_symbols()
        for s, ss in zip(symbols, self._cell.get_symbols()):
            self.assertEqual(s, ss)

    def test_masses(self):
        masses = self._cell.masses
        self._cell.masses = masses
        np.testing.assert_allclose(masses, self._cell.masses)
        self._cell.set_masses(masses)
        self._cell.get_masses()
        np.testing.assert_allclose(masses, self._cell.get_masses())

    def test_magnetic_moments(self):
        magnetic_moments = self._cell.magnetic_moments
        self._cell.magnetic_moments = magnetic_moments
        if magnetic_moments is None:
            self.assertEqual(magnetic_moments, self._cell.magnetic_moments)
        else:
            np.testing.assert_allclose(magnetic_moments, self._cell.magnetic_moments)
        self._cell.get_magnetic_moments()
        self._cell.set_magnetic_moments(magnetic_moments)

    def test_volume(self):
        volume = self._cell.volume
        self._cell.get_volume()
        self.assertTrue(abs(volume - self._cell.get_volume()) < 1e-8)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCell)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
