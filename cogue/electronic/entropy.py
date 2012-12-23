# Entropy of non-interacting electrons
# S=-gk_{\mathrm{B}}\Sigma_i\[f_i \ln f_i + (1-f_i)\ln (1-f_i)\]
# g: degree of degeneracy
# f_i:probability that a state is occpied f_i = \left[1+\exp\left(\frac{E-\mu}{T}\right)\right\]^{-1}

import numpy as np
from cogue.units import Kb

def get_entropy(energies, weights, chemical_potential, temperature):
    mu = chemical_potential
    T = temperature
    S = 0
    for E, w in zip(np.array(energies), weights):
        f = 1.0 / (1 + np.exp((E - mu) / T))
        f = np.extract((f > 1e-10) * (f < 1 - 1e-10), f)
        S += - np.sum(f * np.log(f) + (1 - f) * np.log(1 - f)) * w
    return S

if __name__ == '__main__':
    from cogue.interface.vasp_io import Vasprunxml
    import sys
    
    T = 0.4
    print "Temperature %f K (%f eV)" % (T / Kb, T)
    mu = 8.0724
    vasprun = Vasprunxml(sys.argv[1])
    succeeded = vasprun.parse_eigenvalues()
    eigvals = vasprun.get_eigenvalues()
    kpoints, weights = vasprun.get_kpoints()
    print get_entropy(eigvals, weights, mu, T) * T * 2
