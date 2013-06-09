""" """
import numpy as np
from cogue.crystal.atom import atomic_symbols, atomic_weights

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
