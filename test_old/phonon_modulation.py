"""Test of get_unstable_modulations."""

import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_SETS
from phonopy.interface.vasp import read_vasp

from cogue.task.phonon_relax import get_unstable_modulations

cutoff_eigenvalue = -0.02
supercell_dimension = [2, 2, 2]
cell = read_vasp("POSCAR-unitcell")
phonon = Phonopy(cell, np.diag(supercell_dimension))
force_sets = parse_FORCE_SETS()
phonon.set_displacement_dataset(force_sets)
phonon.produce_force_constants()
phonon.set_mesh(supercell_dimension, is_gamma_center=True)

get_unstable_modulations(
    phonon,
    supercell_dimension,
    cutoff_eigenvalue=cutoff_eigenvalue,
    symmetry_tolerance=0.1,
    max_displacement=0.11,
)
