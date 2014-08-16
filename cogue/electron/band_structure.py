import numpy as np

class BandStructure:
    def __init__(self,
                 paths,
                 cell,
                 eigenvalues,
                 fermi_energy=None):
        self._cell = cell
        self._paths = paths
        self._fermi_energy = fermi_energy
        self._rec_lattice = np.linalg.inv(cell.get_lattice()).T
        self._eigvals = eigenvalues
        self._distances = []
        self._set_distances()

    def write_yaml(self):
        self._write_yaml()

    def _set_distances(self):
        d = 0.0
        for kpoints in self._paths:
            prev_k = kpoints[0]
            dists_path = []
            dists_path.append(d)
            for i, k in enumerate(kpoints[1:]):
                diff = np.array(k) - np.array(prev_k)
                d += np.linalg.norm(np.dot(self._rec_lattice, diff))
                dists_path.append(d)
                prev_k = k
            self._distances.append(dists_path)

    def _write_yaml(self):
        w = open('band.yaml', 'w')
        natom = len(self._cell.get_symbols())
        lattice = self._rec_lattice
        nkpoint = 0
        for kpoints in self._paths:
            nkpoint += len(kpoints)
        w.write("nkpoint: %-7d\n" % nkpoint)
        w.write("npath: %-7d\n" % len(self._paths))
        w.write("natom: %-7d\n" % (natom))
        w.write("reciprocal_lattice:\n")
        for vec, axis in zip(lattice.T, ('a*', 'b*', 'c*')):
            w.write("- [ %12.8f, %12.8f, %12.8f ] # %2s\n" %
                    (tuple(vec) + (axis,)))
        if self._fermi_energy is not None:
            w.write("fermi-energy: %-15.10f\n" % self._fermi_energy)
        w.write("path:\n")
        for i, (kpoints, distances, eigvals) in enumerate(
                zip(self._paths, self._distances, self._eigvals)):
             for j, kpt in enumerate(kpoints):
                w.write("- k-position: [ %12.7f, %12.7f, %12.7f ]\n" %
                        tuple(kpt))
                w.write("  distance: %12.7f\n" % distances[j])
                w.write("  band:\n")
                for k, eig in enumerate(eigvals[j]):
                    w.write("  - # %d\n" % (k + 1))
                    w.write("    eigenvalue: %15.10f\n" % eig)
                w.write("\n")

        
