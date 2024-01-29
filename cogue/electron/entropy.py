"""Entropy related methods."""

# Entropy of non-interacting electrons
# S=-gk_{\mathrm{B}}\Sigma_i\[f_i \ln f_i + (1-f_i)\ln (1-f_i)\]
# g: degree of degeneracy
# f_i:probability that a state is occpied
#     where f_i = \left[1+\exp\left(\frac{E-\mu}{T}\right)\right\]^{-1}

import numpy as np

from cogue.units import Kb


def get_entropy(energies, weights, chemical_potential, temperature):
    """Return entropy."""
    mu = chemical_potential
    T = temperature
    # Degeneracy of electrons
    # Spin components 1 and 2 are stored in tuple as
    # (spin1, spin2) or (spin1,) if no spin polarized.
    # if len(energies) == 1, g = 2 (doubly degenerate)
    g = 3 - len(energies)
    S = 0
    for energies_spin in energies:
        for E, w in zip(np.array(energies_spin), weights):
            f = 1.0 / (1 + np.exp((E - mu) / T))
            f = np.extract((f > 1e-10) * (f < 1 - 1e-10), f)
            S += -np.sum(f * np.log(f) + (1 - f) * np.log(1 - f)) * w
    return S * g


def get_chemical_potential(energies, weights, temperature, num_electrons):
    """Return chemical potential."""
    T = temperature
    emax = np.max(energies)
    emin = np.min(energies)

    for i in range(100):
        mu = (emax + emin) / 2
        n = _get_number_of_electrons(energies, weights, mu, T)
        if abs(n - num_electrons) < 1e-8:
            break
        elif n < num_electrons:
            emin = mu
        else:
            emax = mu

    return mu


def _get_number_of_electrons(energies, weights, chemical_potential, temperature):
    T = temperature
    g = 3 - len(energies)
    mu = chemical_potential
    n = 0
    for energies_spin in energies:
        for E, w in zip(np.array(energies_spin), weights):
            n += np.sum(1.0 / (1 + np.exp((E - mu) / T))) * w
    return n * g


if __name__ == "__main__":
    import sys

    from cogue.interface.vasp_io import Vasprunxml

    T = 0.4
    print("Temperature %f K (%f eV)" % (T / Kb, T))

    vasprun = Vasprunxml(sys.argv[1])
    succeeded = vasprun.parse_eigenvalues()
    eigvals = vasprun.get_eigenvalues()
    kpoints, weights = vasprun.get_kpoints()
    mu = get_chemical_potential(eigvals, weights, T, 16)
    print("Chemical potential: %f" % mu)
    print("Entropy (T*S): %f" % get_entropy(eigvals, weights, mu, T) * T)
