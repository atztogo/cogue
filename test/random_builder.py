import sys
from cogue.calculator.vasp import *
from cogue.crystal.builder import RandomBuilder
from cogue.crystal.cell import get_distance

builder = RandomBuilder(['Si'] * 6 + ['O'] * 12,
                         volume=250)
cell = builder.build()
if not cell:
    print "Cell build failed."
    sys.exit(1)

write_poscar(cell, "POSCAR")

for i, p1 in enumerate(cell.get_points().T):
    for j, p2 in enumerate(cell.get_points().T):
        print "%3d - %3d: %f" % (i+1, j+1,
                                 get_distance(cell.get_lattice(), p1, p2))
