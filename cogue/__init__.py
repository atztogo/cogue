""" """

from cogue.crystal.cell import Cell
from cogue.crystal.symmetry import get_symmetry_dataset
from cogue.controller.autocalc import AutoCalc

def cell(lattice=None,
         points=None,
         symbols=None,
         masses=None,
         numbers=None):
    """ """
    return Cell(lattice=lattice,
                points=points,
                symbols=symbols,
                masses=masses,
                numbers=numbers)

def symmetry(cell, tolerance=1e-5):
    """ """
    return get_symmetry_dataset(cell, tolerance)

def autocalc(name=None, log_name=None, verbose=False):
    """ """
    return AutoCalc(name=name, log_name=log_name, verbose=verbose)

