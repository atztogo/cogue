"""
Spglib interface for cogue
"""

import cogue._spglib as spg
from cogue.crystal.cell import Cell
import numpy as np

def get_symmetry_dataset(cell, tolerance=1e-5):
    """Return crystal symmetry information by analyzing Cell object.

    number: International space group number

    international: International symbol

    hall: Hall symbol

    transformation_matrix:
    Transformation matrix from lattice of input cell to Bravais lattice
    L^bravais = L^original * Tmat

    origin shift: Origin shift in the setting of 'Bravais lattice'

    rotations, translations:
    Rotation matrices and translation vectors
    Space group operations are obtained by
    [(r,t) for r, t in zip(rotations, translations)]

    wyckoffs: Wyckoff letters

    """
    points = cell.get_points().copy()
    lattice = cell.get_lattice().copy()
    numbers = cell.get_numbers()
    keys = ('number',
            'international',
            'hall_number',
            'hall',
            'transformation_matrix',
            'origin_shift',
            'rotations',
            'translations',
            'wyckoffs',
            'equivalent_atoms')
    dataset = {}
    for key, data in zip(keys, spg.get_dataset(lattice, points, numbers, tolerance)):
        dataset[key] = data

    dataset['international'] = dataset['international'].strip()
    dataset['hall'] = dataset['hall'].strip()
    dataset['transformation_matrix'] = np.array(dataset['transformation_matrix'])
    dataset['origin_shift'] = np.array(dataset['origin_shift'])
    dataset['rotations'] = np.array(dataset['rotations'])
    dataset['translations'] = np.array(dataset['translations'])
    letters = "abcdefghijklmnopqrstuvwxyz"
    dataset['wyckoffs'] = [letters[x] for x in dataset['wyckoffs']]
    dataset['equivalent_atoms'] = np.array(dataset['equivalent_atoms'])

    return dataset

def get_crystallographic_cell(cell, tolerance=1e-5):
    points = cell.get_points().copy()
    lattice = cell.get_lattice().copy()
    numbers = cell.get_numbers()
    brv_lattice, brv_points, brv_numbers = \
        spg.get_crystallographic_cell(lattice, points, numbers, tolerance)
    return Cell(lattice=brv_lattice,
                points=brv_points,
                numbers=brv_numbers)

def get_primitive_cell(cell, tolerance=1e-5):
    points = cell.get_points().copy()
    lattice = cell.get_lattice().copy()
    numbers = cell.get_numbers()
    (prim_lattice,
     prim_points,
     prim_numbers) = spg.get_primitive_cell(lattice, points, numbers, tolerance)
    return Cell(lattice=prim_lattice,
                points=prim_points,
                numbers=prim_numbers)
    
