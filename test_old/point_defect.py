import numpy as np

import cogue
import cogue.calculator.vasp as vasp
from cogue.crystal.point_defect import PointDefect
from cogue.crystal.supercell import get_supercell

symbols = ["Cu"] * 4
lattice = [[3.61, 0, 0], [0, 3.61, 0], [0, 0, 3.61]]  # Orthorhombic
points = np.transpose(
    [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
)
cell = cogue.cell(lattice=lattice, points=points, symbols=symbols)
supercell = get_supercell(cell, np.diag([2, 2, 2]))
vacancy_cell = PointDefect(supercell)
vacancy_cell.set_point_vacancy(3)
# To use get_cell restricts the cell to be an object of Cell but not CellBuilder
vasp.write_poscar(vacancy_cell)
