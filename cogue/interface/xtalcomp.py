import cogue._xtalcomp as xcmp
import numpy as np

def compare(cell1,
            cell2,
            tolerance=0.01,
            angle_tolerance=0.1):

    result = xcmp.compare(
        np.array(cell1.lattice.T, dtype='double', order='C'),
        cell1.numbers,
        np.array(cell1.get_points().T, dtype='double', order='C'),
        np.array(cell2.lattice.T, dtype='double', order='C'),
        cell2.numbers,
        np.array(cell2.get_points().T, dtype='double', order='C'),
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

    print(cell1.lattice.T)
    print(cell1.numbers)
    print(cell1.get_points().T)
    print(np.dot(cell1.lattice, cell1.get_points()).T)

    print(cell2.lattice.T)
    print(cell2.numbers)
    print(cell2.get_points().T)
    print(np.dot(cell2.lattice, cell2.get_points()).T)

    print(np.dot(cell1.lattice, cell1.get_points()).T -
          np.dot(cell2.lattice, cell2.get_points()).T)
    for v in (np.dot(cell1.lattice, cell1.get_points()).T -
              np.dot(cell2.lattice, cell2.get_points()).T):
        print(np.dot(v, v) ** (0.5))
    
    print(compare(cell1, cell2, float(sys.argv[3]), float(sys.argv[4])))
