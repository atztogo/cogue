""" """
import numpy as np
from cogue.crystal.atom import atomic_symbols, atomic_weights
import warnings


def symbols2formula(symbols):
    counts = {}
    formula = ""
    for s in set(symbols):
        counts[s] = 0
    for s in symbols:
        counts[s] += 1
    for s in counts:
        formula += "%s%d" % (s, counts[s])
    return formula


def sort_cell_by_symbols(cell):
    symbols = cell.get_symbols()
    compressed_symbols = []
    for s in symbols:
        if not s in compressed_symbols:
            compressed_symbols.append(s)

    atom_order = []
    for cs in compressed_symbols:
        for i, s in enumerate(symbols):
            if s == cs:
                atom_order.append(i)

    if cell.get_magnetic_moments() is not None:
        magmoms = cell.get_magnetic_moments()[atom_order]
    else:
        magmoms = None

    return Cell(lattice=cell.lattice,
                points=(cell.get_points().T)[atom_order].T,
                symbols=[symbols[i] for i in atom_order],
                magmoms=magmoms,
                masses=cell.get_masses()[atom_order])


def get_strained_cells(cell_orig, strains):
    cells = []
    lattice = cell_orig.lattice
    for strain in strains:
        cell = cell_orig.copy()
        if isinstance(strain, int) or isinstance(strain, float):
            cell.set_lattice(lattice * (1 + strain) ** (1.0 / 3))
        else:
            cell.set_lattice(
                np.dot(lattice, np.eye(3) + np.array(strain)))
        cells.append(cell)

    return cells


class Cell(object):
    """ """
    def __init__(self,
                 lattice=None,
                 points=None,
                 symbols=None,
                 magmoms=None,
                 masses=None,
                 numbers=None):

        if lattice is None:
            self._lattice = None
        else:
            self._lattice = np.array(lattice, dtype='double', order='C')

        if points is None:
            self._points = None
        else:
            self._points = np.array(points, dtype='double', order='C')

        if magmoms is None:
            self._magmoms = None
        else:
            self._magmoms = np.array(mogmoms, dtype='double')

        if symbols is None:
            self._symbols = None
        else:
            self._symbols = symbols[:]

        if masses is None:
            self._masses = None
        else:
            self._masses = np.array(masses, dtype='double')

        if numbers is None:
            self._numbers = None
        else:
            self._numbers = np.array(numbers, dtype='intc')

        if self._numbers is None and self._symbols is not None:
            self._set_numbers_from_symbols()

        if self._symbols is None and self._numbers is not None:
            self._set_symbols_from_numbers()

        if self._masses is None:
            self._set_masses_from_numbers()

    @property
    def lattice(self):
        """ """
        return self._lattice.copy()

    @lattice.setter
    def lattice(self, lattice):
        self._lattice = np.array(lattice, dtype='double', order='C')

    def set_lattice(self, lattice):
        """ """
        warnings.warn("set_lattice method is deprecated.", DeprecationWarning)
        self.lattice = lattice

    def get_lattice(self):
        """ """
        warnings.warn("get_lattice method is deprecated.", DeprecationWarning)
        return self.lattice

    @property
    def volume(self):
        """ """
        return np.linalg.det(self._lattice)

    def get_volume(self):
        """ """
        warnings.warn("get_volume method is deprecated.", DeprecationWarning)
        return self.volume

    @property
    def points(self):
        """ """
        return self._points.copy()

    @points.setter
    def points(self, points):
        """ """
        self._points = np.array(points, dtype='double', order='C')

    def set_points(self, points):
        """ """
        warnings.warn("set_points method is deprecated.", DeprecationWarning)
        self.points = points

    def get_points(self):
        """ """
        warnings.warn("get_points method is deprecated.", DeprecationWarning)
        return self.points

    @property
    def symbols(self):
        """ """
        return self._symbols[:]

    @symbols.setter
    def symbols(self, symbols):
        """ """
        self._symbols = symbols[:]
        self._set_numbers_from_symbols()
        self._set_masses_from_numbers()

    def set_symbols(self, symbols):
        """ """
        warnings.warn("set_symbols method is deprecated.", DeprecationWarning)
        self.symbols = symbols

    def get_symbols(self):
        """ """
        warnings.warn("get_symbols method is deprecated.", DeprecationWarning)
        return self.symbols

    @property
    def masses(self):
        """ """
        return self._masses.copy()

    @masses.setter
    def masses(self, masses):
        """ """
        self._masses = np.array(masses, dtype='double')

    def set_masses(self, masses):
        """ """
        warnings.warn("set_masses method is deprecated.", DeprecationWarning)
        self.masses = masses

    def get_masses(self):
        """ """
        warnings.warn("get_masses method is deprecated.", DeprecationWarning)
        return self.masses

    @property
    def magnetic_moments(self):
        """ """
        if self._magmoms is None:
            return None
        else:
            return self._magmoms.copy()

    @magnetic_moments.setter
    def magnetic_moments(self, magmoms):
        """ """
        if magmoms is None:
            self._magmoms = None
        else:
            self._magmoms = np.array(magmoms, dtype='double')

    def set_magnetic_moments(self, magmoms):
        """ """
        warnings.warn("set_magnetic_moments method is deprecated.",
                      DeprecationWarning)
        self.magnetic_moments = magmoms

    def get_magnetic_moments(self):
        """ """
        warnings.warn("get_magnetic_moments method is deprecated.",
                      DeprecationWarning)
        return self.magnetic_moments

    @property
    def numbers(self):
        return self._numbers.copy()

    @numbers.setter
    def numbers(self, numbers):
        """ """
        self._numbers = np.array(numbers, dtype='intc')
        self._set_symbols_from_numbers()
        self._set_masses_from_numbers()

    def set_numbers(self, numbers):
        """ """
        warnings.warn("set_numbers method is deprecated.", DeprecationWarning)
        self.numbers = numbers

    def get_numbers(self):
        """ """
        warnings.warn("get_numbers method is deprecated.", DeprecationWarning)
        return self.numbers

    def copy(self):
        """ """
        return Cell(lattice=self._lattice,
                    points=self._points,
                    numbers=self._numbers,
                    magmoms=self._magmoms,
                    masses=self._masses)

    def get_yaml_lines(self):
        lines = []
        lines.append("lattice:")
        for v, a in zip(self._lattice.T, ('a', 'b', 'c')):
            lines.append("- [ %22.16f, %22.16f, %22.16f ] # %s" %
                         (v[0], v[1], v[2], a))

        lines.append("points:")
        for i, (s, v, m) in enumerate(zip(self._symbols,
                                          self._points.T,
                                          self._masses)):
            lines.append("- symbol: %-2s # %d" % (s, i + 1))
            lines.append("  coordinates: [ %19.16f, %19.16f, %19.16f ]" %
                         tuple(v))
            lines.append("  mass: %f" % m)

        return lines

    def __str__(self):
        return "\n".join(self.get_yaml_lines())

    def _set_numbers_from_symbols(self):
        self._numbers = np.array([atomic_symbols[s] for s in self._symbols],
                                 dtype='intc')

    def _set_symbols_from_numbers(self):
        self._symbols = [atomic_weights[x][0] for x in self._numbers]

    def _set_masses_from_numbers(self):
        self._masses = np.array([atomic_weights[x][3] for x in self._numbers],
                                dtype='double')
