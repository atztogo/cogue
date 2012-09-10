import cogue._xtalcomp as xcmp
import numpy as np

def compare(cell1,
            cell2,
            tolerance=0.01,
            angle_tolerance=0.1):

    result = xcmp.compare(cell1.get_lattice().T.copy(),
                          cell1.get_numbers(),
                          cell1.get_points().T.copy(),
                          cell2.get_lattice().T.copy(),
                          cell2.get_numbers(),
                          cell2.get_points().T.copy(),
                          tolerance,
                          angle_tolerance)
    if result == 0:
        return False
    else:
        return True

if __name__ == '__main__':
    import sys
    from cogue.calculator.vasp import read_poscar

    cell1 = read_poscar(sys.argv[1])
    cell2 = read_poscar(sys.argv[2])

    print cell1.get_lattice().T.copy()
    print cell1.get_numbers()
    print cell1.get_points().T.copy()

    print cell2.get_lattice().T.copy()
    print cell2.get_numbers()
    print cell2.get_points().T.copy()
    
    print compare(cell1, cell2, float(sys.argv[3]))
