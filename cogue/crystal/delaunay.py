import numpy as np

#
# Delaunay reduction (originally from phonopy)
#    
def get_Delaunay_reduction(lattice, tolerance):
    extended_bases = np.zeros((4, 3), dtype=float)
    extended_bases[:3,:] = np.transpose(lattice)
    extended_bases[3]  = -np.sum(lattice, axis=1)

    for i in range(100):
        if reduce_bases(extended_bases, tolerance):
            break
    if i == 99:
        print("Delaunary reduction failed.")

    shortest3 = get_shortest_bases(extended_bases, tolerance)

    return shortest3.T.copy()

def reduce_bases(extended_bases, tolerance):
    metric = np.dot(extended_bases, extended_bases.T)
    for i in range(4):
        for j in range(i+1, 4):
            if metric[i, j] > tolerance:
                for k in range(4):
                    if (not k == i) and (not k == j):
                        extended_bases[k] += extended_bases[i]
                extended_bases[i] = -extended_bases[i]
                extended_bases[j] =  extended_bases[j]
                return False

    # All non diagonal elements of metric tensor is negative.
    # Reduction is completed.
    return True

def get_shortest_bases(extended_bases, tolerance):

    def mycmp(x, y):
        return cmp(np.vdot(x, x), np.vdot(y, y))

    basis = np.zeros((7, 3), dtype=float)
    basis[:4] = extended_bases
    basis[4]  = extended_bases[0] + extended_bases[1]
    basis[5]  = extended_bases[1] + extended_bases[2]
    basis[6]  = extended_bases[2] + extended_bases[0]
    # Sort bases by the lengthes (shorter is earlier)
    basis = sorted(basis, cmp=mycmp)
    
    # Choose shortest and linearly independent three bases
    # This algorithm may not be perfect.
    for i in range(7):
        for j in range(i+1, 7):
            for k in range(j+1, 7):
                if abs(np.linalg.det([basis[i],
                                      basis[j],
                                      basis[k]])) > tolerance:
                    return np.array([basis[i],
                                     basis[j],
                                     basis[k]])

    print("Delaunary reduction is failed.")
    return basis[:3]
