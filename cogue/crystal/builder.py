import sys
import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.atom import atomic_symbols, atomic_weights

class CellBuilder:
    def __init__(self, cell):
        self._points = cell.get_points()
        self._symbols = cell.get_symbols()
        self._magmoms = cell.get_magnetic_moments()
        self._masses = cell.get_masses()
        self._numbers = cell.get_numbers()
        self._lattice = cell.get_lattice()
        
    def push(self,
             point=None,
             symbol=None,
             magmom=None,
             mass=None,
             number=None):
        if self._magmoms is not None and magmom is None:
            sys.stderr.write("Magmoms has to be set.\n")
        elif point is not None and (symbol or number):
            self._points = np.append(self._points, [[point[0]],
                                                    [point[1]],
                                                    [point[2]]], axis=1)
            if symbol:
                self._symbols.append(symbol)
            if number:
                self._numbers = np.append(self._numbers, number)
            if not symbol:
                self._symbols.append(atomic_weights[number][0])
            if not number:
                self._numbers = np.append(self._numbers, atomic_symbols[symbol])
            if atomic_symbols[self._symbols[-1]] != self._numbers[-1]:
                sys.stderr.write(
                    "Symbol and number don't match. Symbol is taken.\n")
                self._numbers[-1] = atomic_symbols[symbol]
            if mass:
                self._masses = np.append(self._masses, mass)
            else:
                self._masses = np.append(self._masses,
                                         atomic_weights[self._numbers[-1]][3])

            if magmom is not None:
                if self._magmoms is None:
                    sys.stderr.write(
                        "Magmoms is not defined.\n")
                else:
                    self._points = np.append(self._magmoms, [magmom], axis=0)
        else:
            sys.stderr.write(
                "At least a pair of point and symbol or "
                "a pair of point and number have to be set.\n")

    def pop(self, index=None):
        if index is None:
            i = len(self._symbols) - 1
        else:
            i = index
        self._points = np.delete(self._points, i, axis=1)
        self._symbols.pop(i)
        if self._magmoms is not None:
            self._magmoms = np.delete(self._magmoms, i)
        self._masses = np.delete(self._masses, i)
        self._numbers = np.delete(self._numbers, i)

    def get_cell(self):
        return Cell(lattice=self._lattice,
                    magmoms=self._magmoms,
                    masses=self._masses,
                    numbers=self._numbers,
                    points=self._points)
                    

