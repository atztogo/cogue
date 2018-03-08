import sys
import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.symmetry import (get_symmetry_dataset,
                                    get_crystallographic_cell)

##############################
# Generally usable functions #
##############################
def reduce_points(tmat, cell, tolerance=1e-5):
    lattice = cell.lattice
    points_prim = []
    symbols_prim = []
    masses_prim = []
    points = cell.get_points()
    symbols = cell.get_symbols()
    masses = cell.get_masses()
    magmoms = cell.get_magnetic_moments()
    if magmoms==None:
        magmoms_prim = None
    else:
        magmoms_prim = []

    for i, p in enumerate(np.dot(np.linalg.inv(tmat), points).T):
        is_different = True
        for p_prim in points_prim:
            diff = p_prim - p
            diff -= diff.round()
            if (np.linalg.norm(np.dot(lattice, diff)) < tolerance).all():
                is_different = False
                break
        if is_different:
            points_prim.append(p - np.floor(p))
            symbols_prim.append(symbols[i])
            masses_prim.append(masses[i])
            if magmoms != None:
                magmoms_prim.append(magmoms[i])

    return Cell(lattice=np.dot(cell.lattice, tmat),
                points=np.transpose(points_prim),
                symbols=symbols_prim,
                magmoms=magmoms_prim,
                masses=masses_prim)

def get_primitive(cell, tolerance=1e-5):
    # spglib returns R-centred lattice for Rhombohedrals
    std_cell = get_crystallographic_cell(cell, tolerance)
    sym_dataset = get_symmetry_dataset(std_cell)
    spg_symbol = sym_dataset['international'][0]
    if spg_symbol == 'F':
        std_cell = _fc2prim(std_cell)
    elif spg_symbol == 'I':
        std_cell = _bc2prim(std_cell)
    elif spg_symbol == 'A':
        std_cell = _abc2prim(std_cell)
    elif spg_symbol == 'B':
        std_cell = _bbc2prim(std_cell)
    elif spg_symbol == 'C':
        std_cell = _cbc2prim(std_cell)
    return std_cell

def _fc2prim(cell):
    tmat = [[ 0.0, 0.5, 0.5],
            [ 0.5, 0.0, 0.5],
            [ 0.5, 0.5, 0.0]]
    return reduce_points(tmat, cell)

def _bc2prim(cell):
    tmat = [[-0.5, 0.5, 0.5],
            [ 0.5,-0.5, 0.5],
            [ 0.5, 0.5,-0.5]]
    return reduce_points(tmat, cell)

def _cbc2prim(cell):
    tmat = [[ 0.5, 0.5, 0.0],
            [-0.5, 0.5, 0.0],
            [ 0.0, 0.0, 1.0]]
    return reduce_points(tmat, cell)

def _abc2prim(cell):
    tmat = [[ 0.0, 0.5, 0.5],
            [ 0.0,-0.5, 0.5],
            [ 1.0, 0.0, 0.0]]
    return reduce_points(tmat, cell)

def _bbc2prim(cell):
    tmat = [[ 0.5, 0.0, 0.5],
            [ 0.5, 0.0,-0.5],
            [ 0.0, 1.0, 0.0]]
    return reduce_points(tmat, cell)

#########################
# Phonopy Atoms to Cell #
#########################
def atoms2cell(phonopy_cell):
    return Cell(lattice=phonopy_cell.get_cell().T,
                points=phonopy_cell.get_scaled_positions().T,
                masses=phonopy_cell.get_masses(),
                magmoms=phonopy_cell.get_magnetic_moments(),
                symbols=phonopy_cell.symbols)

#########################
# Cell to Phonopy Atoms #
#########################
def cell2atoms(cell):
    try:
        from phonopy.structure.atoms import PhonopyAtoms as Atoms
    except ImportError:
        print("You need to install phonopy.")
        sys.exit(1)
    return Atoms(cell=cell.lattice.T,
                 scaled_positions=cell.get_points().T,
                 masses=cell.get_masses(),
                 magmoms=cell.get_magnetic_moments(),
                 symbols=cell.get_symbols())

#######################
# Writers and readers #
#######################

#
# yaml
#
def read_yaml(filename):
    import yaml
    data = yaml.load(open(filename))

    if 'status' in data:
        if (not data['status'] == 'terminate' and
            'lattice' in data and
            'points' in data and
            'symbols' in data):
            lattice = np.transpose(data['lattice'])
            points = np.transpose(data['points'])
            symbols = data['symbols']
            return Cell(lattice=lattice,
                        points=points,
                        symbols=symbols)
    return None



if __name__ == '__main__':
    pass
