import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.converter import reduce_points
from cogue.crystal.utility import get_lattice_parameters
from cogue.crystal.symmetry import get_symmetry_dataset

def get_supercell(cell, supercell_matrix, tolerance=1e-5):
    """
    Build supercell with supercell matrix as follows:
    1) Create supercell with surrounding simple lattice.
    2) Trim the surrounding supercell with the target lattice.

    """
    
    smat = np.array(supercell_matrix)
    frame = get_smallest_surrounding_lattice_multiplicities(smat)
    surrounding_cell = get_simple_supercell(frame, cell)
    if (np.diag(np.diagonal(smat)) == smat).all():
        return surrounding_cell

    transformation_matrix = np.array([smat[0] / float(frame[0]),
                                      smat[1] / float(frame[1]),
                                      smat[2] / float(frame[2])])
    supercell = reduce_points(transformation_matrix,
                              surrounding_cell,
                              tolerance)
    return supercell

def get_simple_supercell(multi, cell):
    lattice = cell.get_lattice()
    points = cell.get_points()
    masses = cell.get_masses()
    magmoms = cell.get_magnetic_moments()
    symbols = cell.get_symbols()
    
    points_scell = []
    symbols_scell = []
    masses_scell = []
    if magmoms == None:
        magmoms_scell = None
    else:
        magmoms_scell = []
    for l, pos in enumerate(points.T):
        for i in range(multi[2]):
            for j in range(multi[1]):
                for k in range(multi[0]):
                    points_scell.append([(pos[0] + k) / multi[0],
                                         (pos[1] + j) / multi[1],
                                         (pos[2] + i) / multi[2]])
                    symbols_scell.append(symbols[l])
                    masses_scell.append(masses[l])
                    if not magmoms == None:
                        magmoms_scell.append(magmoms[l])

    return Cell(lattice=np.dot(lattice, np.diag(multi)),
                points=np.transpose(points_scell),
                symbols=symbols_scell,
                masses=masses_scell,
                magmoms=magmoms_scell)

def get_smallest_surrounding_lattice_multiplicities(supercell_matrix):
    """Build a frame surrounding supercell lattice"""

    m = np.array(supercell_matrix)
    axes = np.array([[ 0, 0, 0],
                     m[:,0],
                     m[:,1],
                     m[:,2],
                     m[:,1] + m[:,2],
                     m[:,2] + m[:,0],
                     m[:,0] + m[:,1],
                     m[:,0] + m[:,1] + m[:,2]])
    frame = [max(axes[:, i]) - min(axes[:, i]) for i in (0, 1, 2)]
    return frame

def estimate_supercell_matrix(cell,
                              max_num_atoms=120,
                              symprec=1e-5):
    """Estimate supercell matrix from conventional cell

    Supercell matrix is estimated from basis vector lengths and
    maximum number of atoms accepted. The input cell must be already
    standarized. For triclinic, monoclinic, and orthorhombic cells,
    basis vectors of a, b, c are freely multiplied, but for tetragonal
    and hexagonal cells, multiplicities of a and b are the same, and
    for cubic cell, those of a, b, c are to be found the same.  The
    supercell matrix is always returned as a diagonal matrix.

    """

    dataset = get_symmetry_dataset(cell)
    spg_num = dataset['number']
    num_atoms = len(cell.get_numbers())
    lengths_orig = get_lattice_parameters(cell.get_lattice())
    lengths = get_lattice_parameters(dataset['std_lattice'])

    assert (np.abs(lengths_orig - lengths) < symprec).all()

    if spg_num <= 74: # Triclinic, monoclinic, and orthorhombic
        multi = _get_multiplicity_abc(num_atoms, lengths, max_num_atoms)
    elif spg_num <= 194: # Tetragonal and hexagonal
        multi = _get_multiplicity_ac(num_atoms, lengths, max_num_atoms)
    else: # Cubic
        multi = _get_multiplicity_a(num_atoms, lengths, max_num_atoms)

    smat = np.eye(3, dtype='intc')
    for i in range(3):
        smat[i, i] = multi[i]

    return smat

def _get_multiplicity_abc(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = [1, 1, 1]

    for i in range(max_iter):
        l_super = np.multiply(multi, lengths)
        min_index = np.argmin(l_super)
        multi[min_index] += 1
        if num_atoms * np.prod(multi) > max_num_atoms:
            multi[min_index] -= 1

    return multi

def _get_multiplicity_ac(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = [1, 1]
    a = lengths[0]
    c = lengths[2]

    for i in range(max_iter):
        l_super = np.multiply(multi, [a, c])
        min_index = np.argmin(l_super)
        multi[min_index] += 1
        if num_atoms * multi[0] ** 2 * multi[1] > max_num_atoms:
            multi[min_index] -= 1

    return [multi[0], multi[0], multi[1]]

def _get_multiplicity_a(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = 1
    a = lengths[0]

    for i in range(max_iter):
        l_super = multi * a
        multi += 1
        if num_atoms * multi ** 3 > max_num_atoms:
            multi -= 1

    return [multi, multi, multi]

