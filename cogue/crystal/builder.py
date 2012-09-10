import random
import numpy as np
from cogue.crystal.cell import Cell, get_distance

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


if __name__ == '__main__':
    import sys
    from cogue.calculator.vasp import *

    builder = RandomBuilder(['Si'] * 6 + ['O'] * 12,
                             volume=250)
    cell = builder.build()
    if not cell:
        print "Cell build failed."
        sys.exit(1)

    write_poscar(cell)

    for i, p1 in enumerate(cell.get_points().T):
        for j, p2 in enumerate(cell.get_points().T):
            print "%3d - %3d: %f" % (i+1, j+1,
                                     get_distance(cell.get_lattice(), p1, p2))
