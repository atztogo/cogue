import numpy as np

from cogue.crystal.delaunay import get_Delaunay_reduction


def get_pair_distances(cell, tolerance=1e-5):
    positions = cell.get_points().T
    num_atom = len(positions)
    distances = np.zeros((num_atom, num_atom), dtype=float)
    for i, p1 in enumerate(positions):
        for j, p2 in enumerate(positions):
            distances[i, j] = get_distance(cell.get_lattice(), p1, p2)
    return distances


def get_distance(lattice, p1, p2, tolerance=1e-5):
    """Return shortest distance between a pair of atoms in PBC"""
    shortest_bases = get_Delaunay_reduction(lattice, tolerance)
    diff = np.dot(
        np.dot(np.linalg.inv(shortest_bases), lattice), np.array(p1) - np.array(p2)
    )
    diff -= diff.round()

    distances = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                distances.append(
                    np.linalg.norm(np.dot(shortest_bases, diff + np.array([i, j, k])))
                )
    return min(distances)
