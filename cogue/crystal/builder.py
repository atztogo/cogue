import random
import sys
import numpy as np
from cogue.crystal.cell import Cell, get_distance
from cogue.crystal.atom import atomic_symbols, atomic_weights

class CellBuilder(Cell):
    def __init__(self, cell):
        Cell.__init__(self,
                      lattice=cell.get_lattice(),
                      points=cell.get_points(),
                      numbers=cell.get_numbers(),
                      magmoms=cell.get_magnetic_moments(),
                      masses=cell.get_masses())
        
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
        return self.copy()

class RandomBuilder:
    def __init__(self,
                 symbols,
                 volume=None,
                 min_distance=None,
                 max_distance=None):

        self._symbols = symbols
        self._min_distance = min_distance
        if not volume:
            self._volume = 25.0 * len(symbols)
        else:
            self._volume = volume

        if not max_distance:
            self._max_distance = self._volume ** (1.0 / 3) * 3 ** (1.0 / 2) / 2
        else:
            self._max_distance = max_distance

        if not min_distance:
            self._min_distance = (self._volume / len(self._symbols) / (4.0 / 3 * np.pi)) ** (1.0 / 3) * 1.3
        else:
            self._min_distance = min_distance

    def build(self):
        for i in range(100):
            cell = self._random()
            if cell:
                break
            else:
                self._min_distance *= 0.99

        return self._shuffle(cell)

    def get_min_distance(self):
        return self._min_distance

    def set_min_distance(self, min_distance):
        self._min_distance = min_distance

    def get_max_distance(self):
        return self._max_distance

    def set_max_distance(self, max_distance):
        self._max_distance = max_distance

    def _shuffle(self, cell):
        indices = range(len(cell.get_symbols()))
        random.shuffle(indices)
        points = np.zeros(cell.get_points().shape, dtype=float)
        for i, j in enumerate(indices):
            points[:,i] = cell.get_points()[:,j]
            
        return Cell(lattice = cell.get_lattice(),
                    points = points,
                    symbols = self._symbols)

    def _random(self):
        lattice = np.eye(3, dtype=float) * self._volume ** (1.0/3)
        points = [[0., 0., 0.]]
        for i in range(len(self._symbols) - 1):
            x = self._get_random_point(lattice, points)
            if x:
                points.append(x)
            if len(points) == len(self._symbols):
                break

        if len(points) < len(self._symbols):
            return None
        else:
            return Cell(lattice = lattice,
                        points = np.transpose(points),
                        symbols = self._symbols)

    def _get_random_point(self, lattice, points):
        attempt = 0
        while True:
            x = [random.random(),
                 random.random(),
                 random.random()]

            all_ok = True
            for p in points:
                d = get_distance(lattice, x, p)
                if (d < self._min_distance or
                    d > self._max_distance):
                    all_ok = False
                    break

            if all_ok:
                return x

            attempt += 1
            if attempt == 20:
                return False


