import numpy as np

import cogue
import cogue.calculator.vasp as vasp
from cogue.crystal.builder import CellBuilder

symbols = ["Si"] * 2 + ["O"] * 4
lattice = [[4.65, 0, 0], [0, 4.75, 0], [0, 0, 3.25]]  # Orthorhombic
points = np.transpose(
    [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
)
cell = cogue.cell(lattice=lattice, points=points, symbols=symbols)

builder = CellBuilder(cell)
builder.push(point=[0.1, 0.1, 0.1], symbol="Si")
builder.pop(4)
# To use get_cell restricts the cell to be an object of Cell but not CellBuilder
vasp.write_poscar(builder.get_cell())
