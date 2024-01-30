"""Random structure builder."""

import sys

from cogue.crystal.pair_distance import get_distance
from cogue.crystal.random_builder import RandomBuilder
from cogue.interface.vasp_io import write_poscar

builder = RandomBuilder(["Si"] * 6 + ["O"] * 12, volume=250)
cell = builder.build()
if not cell:
    print("Cell build failed.")
    sys.exit(1)

write_poscar(cell, "POSCAR")

for i, p1 in enumerate(cell.get_points().T):
    for j, p2 in enumerate(cell.get_points().T):
        print(
            "%3d - %3d: %f" % (i + 1, j + 1, get_distance(cell.get_lattice(), p1, p2))
        )
