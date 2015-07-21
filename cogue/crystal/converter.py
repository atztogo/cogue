import sys
import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.symmetry import get_symmetry_dataset, get_crystallographic_cell

##############################
# Generally usable functions #
##############################
def reduce_points(tmat, cell, tolerance=1e-5):
    lattice = cell.get_lattice()
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

    return Cell(lattice=np.dot(cell.get_lattice(), tmat),
                points=np.transpose(points_prim),
                symbols=symbols_prim,
                magmoms=magmoms_prim,
                masses=masses_prim)

#####################
# Utility functions #
#####################
def frac2val(string):
    if '/' in string:
        num, denom = [float(x) for x in string.split('/')]
        return num / denom
    else:
        return float(string)

def get_angles(lattice):
    a, b, c = get_lattice_parameters(lattice)
    alpha = np.arccos(np.vdot(lattice[:,1], lattice[:,2]) / b / c)
    beta  = np.arccos(np.vdot(lattice[:,2], lattice[:,0]) / c / a)
    gamma = np.arccos(np.vdot(lattice[:,0], lattice[:,1]) / a / b)
    return alpha / np.pi * 180, beta / np.pi * 180, gamma / np.pi * 180

def lattice2cartesian(a, b, c, alpha, beta, gamma):
    """
    The conversion refers the wikipedia,
    http://en.wikipedia.org/wiki/Fractional_coordinates
    """
    cg = np.cos(gamma / 180 * np.pi)
    cb = np.cos(beta / 180 * np.pi)
    ca = np.cos(alpha / 180 * np.pi)
    sg = np.sin(gamma / 180 * np.pi)
    L = np.zeros((3, 3))
    L[0, 0] = a
    L[0, 1] = b * cg
    L[0, 2] = c * cb
    L[1, 1] = b * sg
    L[1, 2] = c * (ca - cb * cg) / sg
    L[2, 2] = c * np.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2 +
                          2 * ca * cb * cg) / sg
    return L

def get_lattice_parameters(lattice):
    return np.sqrt(np.dot(lattice.T, lattice).diagonal())

def get_oriented_lattice(lattice):
    a, b, c = get_lattice_parameters(lattice)
    alpha, beta, gamma = get_angles(lattice)
    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    a1 = a
    a2 = 0.0
    a3 = 0.0
    b1 = np.cos(gamma)
    b2 = np.sin(gamma)
    b3 = 0.0
    c1 = np.cos(beta)
    c2 = (2 * np.cos(alpha) + b1**2 + b2**2 - 2 * b1 * c1 - 1) / (2 * b2)
    c3 = np.sqrt(1 - c1**2 - c2**2)
    lattice = np.zeros((3, 3), dtype=float)
    lattice[0, 0] = a
    lattice[:,1] = np.array([b1, b2, b3]) * b
    lattice[:,2] = np.array([c1, c2, c3]) * c
    return lattice

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
        from phonopy.structure.atoms import Atoms
    except ImportError:
        print "You need to install phonopy."
        exit(1)
    return Atoms(cell=cell.get_lattice().T,
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
            return Cell(lattice = lattice,
                        points = points,
                        symbols = symbols)
    return None

            
                                                                        
if __name__ == '__main__':
    pass
