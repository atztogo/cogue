""" """
import numpy as np
from cogue.crystal.delaunay import get_Delaunay_reduction
from cogue.crystal.atom import atomic_symbols, atomic_weights

def get_pair_distances(cell, tolerance=1e-5):
    positions = cell.get_points().T
    num_atom = len(positions)
    distances = np.zeros((num_atom, num_atom), dtype=float)
    for i, p1 in enumerate(positions):
        for j, p2 in enumerate(positions):
            distances[i, j] = get_distance(cell.get_lattice(), p1, p2)
    return distances

def get_distance(lattice, p1, p2, tolerance=1e-5):
    """Return shortest distance between a pair of atoms in PBC
    """
    shortest_bases = get_Delaunay_reduction(lattice, tolerance)
    diff = np.dot(np.dot(np.linalg.inv(shortest_bases), lattice),
                  np.array(p1) - np.array(p2))
    diff -= diff.round()

    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append(np.linalg.norm(
                        np.dot(shortest_bases,
                               diff + np.array([i, j, k]))))
    return min(distances)


class Cell:
    """ """
    def __init__(self,
                 lattice=None,
                 points=None,
                 symbols=None,
                 magmoms=None,
                 masses=None,
                 numbers=None):

        if lattice == None:
            self._lattice = None
        else:
            self._lattice = np.array(lattice, dtype=float).copy()
            
        if points == None:
            self._points = None
        else:
            self._points = np.array(points, dtype=float).copy()

        if magmoms == None:
            self._magmoms = None
        else:
            self._magmoms = np.array(mogmoms, dtype=float).copy()

        if not symbols:
            self._symbols = None
        else:
            self._symbols = symbols[:]

        if masses == None:
            self._masses = None
        else:
            self._masses = np.array(masses, dtype=float).copy()

        if numbers == None:
            self._numbers = None
        else:
            self._numbers = np.array(numbers).copy()

        if self._numbers == None and self._symbols:
            self._set_numbers_from_symbols()
            
        if not self._symbols and not self._numbers == None:
            self._set_symbols_from_numbers()

        if self._masses == None:
            self._set_masses_from_numbers()

    def _set_numbers_from_symbols(self):
        self._numbers = np.array(
            [atomic_symbols[s] for s in self._symbols])

    def _set_symbols_from_numbers(self):
        self._symbols = [atomic_weights[x][0] for x in self._numbers]

    def _set_masses_from_numbers(self):
        self._masses = np.array(
            [atomic_weights[x][3] for x in self._numbers])

    def set_lattice(self, lattice):
        """ """
        self._lattice = np.array(lattice, dtype=float).copy()

    def get_lattice(self):
        """ """
        return self._lattice

    def get_volume(self):
        """ """
        return np.linalg.det(self._lattice)

    def set_points(self, points):
        """ """
        self._points = np.array(points, dtype=float).copy()

    def get_points(self):
        """ """
        return self._points

    def set_symbols(self, symbols):
        """ """
        self._symbols = symbols[:]
        self._set_numbers_from_symbols()
        self._set_masses_from_numbers()
        
    def get_symbols(self):
        """ """
        return self._symbols

    def set_masses(self, masses):
        """ """
        self._masses = np.array(masses, dtype=float).copy()

    def get_masses(self):
        """ """
        return self._masses

    def set_magnetic_moments(self, magmoms):
        """ """
        self._magmoms = np.array(magmoms, dtype=float).copy()

    def get_magnetic_moments(self):
        """ """
        return self._magmoms

    def set_numbers(self, numbers):
        """ """
        self._numbers = np.array(numbers).copy()
        self._set_symbols_from_numbers()
        self._set_masses_from_numbers()

    def get_numbers(self):
        """ """
        return self._numbers

    def copy(self):
        """ """
        return Cell(lattice=self._lattice,
                    points=self._points,
                    numbers=self._numbers,
                    magmoms=self._magmoms,
                    masses=self._masses)
