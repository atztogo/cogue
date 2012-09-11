"""
Symmetry functions
"""

import cogue.crystal.spglib as spglib

def get_symmetry_dataset(cell, tolerance=1e-5):
    return spglib.get_symmetry_dataset(cell, tolerance)

def get_crystallographic_cell(cell, tolerance=1e-5):
    return spglib.get_crystallographic_cell(cell, tolerance)

def get_primitive_cell(cell, tolerance=1e-5):
    return spglib.get_primitive_cell(cell, tolerance)
