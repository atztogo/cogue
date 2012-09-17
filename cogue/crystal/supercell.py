import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.converter import reduce_points

def get_supercell(cell, supercell_matrix, tolerance=1e-5):
    """
    Build supercell with supercell matrix as follows:
    1) Create supercell with surrounding simple lattice.
    2) Trim the surrounding supercell with the target lattice.

    """
    
    smat = np.array(supercell_matrix)
    frame = get_smallest_surrounding_lattice_multiplicities(smat)
    surrounding_cell = get_simple_supercell(frame, cell)
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
