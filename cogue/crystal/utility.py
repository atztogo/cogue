import numpy as np
from cogue.crystal.symmetry import (get_symmetry_dataset,
                                    get_crystallographic_cell)

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

def klength2mesh(k_length, lattice):
    """Convert length to mesh in k-point sampling

    This conversion follows VASP manual.

    """
    rec_lattice = np.linalg.inv(lattice).T
    rec_lat_lengths = np.sqrt(np.diagonal(np.dot(rec_lattice.T, rec_lattice)))
    k_mesh = (rec_lat_lengths * k_length + 0.5).astype(int)
    return np.maximum(k_mesh, [1, 1, 1])

def get_Z(numbers):
    count = {}
    for n in numbers:
        if n in count:
            count[n] += 1
        else:
            count[n] = 1
    values  = list(count.values())

    try:
        from math import gcd
    except ImportError:
        from fractions import gcd

    x = values[0]
    for v in values[1:]:
        x = gcd(x, v)

    return x
