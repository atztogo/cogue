from cogue.crystal.symmetry import \
    get_symmetry_dataset, get_crystallographic_cell
from cogue.crystal.converter import atoms2cell
from cogue.crystal.cell import Cell
from cogue.crystal.xtalcomp import compare as xtal_compare
import numpy as np

class PhononModulation:
    def __init__(self,
                 phonon,
                 qpoint,
                 band_indices,
                 modulation_dimension,
                 ndiv=180,
                 symmetry_tolerance=0.05,
                 max_displacement=0.1):
        self._phonon = phonon
        self._qpoint = qpoint
        self._band_indices = band_indices
        self._modulation_dimension = modulation_dimension
        self._ndiv = ndiv
        self._symmetry_tolerance = symmetry_tolerance
        self._max_displacement = max_displacement
        self._vectors = []
        self._points_on_sphere = []

        self._run()
        
    def get_modulation_cells(self):
        return [self._get_cell_with_modulation(m)
                for m in self.get_modulations()]

    def get_modulations(self):
        return [self._get_modulation(pt) / ph * a
                for (pt, ph, a) in self._points_on_sphere]

    def get_points_on_sphere(self):
        return self._points_on_sphere

    def get_vectors(self):
        return self._vectors

    def get_supercell(self):
        return self._supercell

    def _run(self):
        self._set_vectors_and_supercell()
        max_num_op = 0
        best_cells = []
        best_spacegroup_types = []
        points_on_sphere = []
        phase_shifts = self._get_phase_shifts_at_lattice_points()
        for point in self._get_phases():
            modulation = self._get_modulation(point)

            for phase in phase_shifts:
                amplitude = self._get_normalize_amplitude(modulation / phase)
                modcell = self._get_cell_with_modulation(
                    modulation / phase * amplitude)
                symmetry = get_symmetry_dataset(
                    modcell, tolerance=self._symmetry_tolerance)
    
                num_op = len(symmetry['rotations'])
                if num_op > max_num_op:
                    max_num_op = num_op
                    best_cells = [modcell]
                    best_spacegroup_types = [symmetry['number']]
                    points_on_sphere = [[point, phase, amplitude]]

                elif num_op == max_num_op:
                    if symmetry['number'] in best_spacegroup_types:
                        cell_in_best_cells = False
                        for bc in best_cells:
                            if xtal_compare(
                                bc,
                                modcell,
                                tolerance=self._symmetry_tolerance,
                                angle_tolerance=1.0):
                                cell_in_best_cells = True
                                break
                                
                        if not cell_in_best_cells:
                            best_cells.append(modcell)
                            points_on_sphere.append([point, phase, amplitude])
                    else:
                        best_cells.append(modcell)
                        points_on_sphere.append([point, phase, amplitude])
                        best_spacegroup_types.append(symmetry['number'])

        self._points_on_sphere = points_on_sphere

    def _set_vectors_and_supercell(self):
        phonon_modes = [[self._qpoint, i, 1, 0]
                        for i in self._band_indices]
        self._phonon.set_modulations(self._modulation_dimension,
                                     phonon_modes)
        modulations, supercell = self._phonon.get_delta_modulations()
        self._vectors = [delta.T for delta in modulations]
        self._supercell = atoms2cell(supercell)
        self._lattice = self._supercell.get_lattice()
        self._lattice_inv = np.linalg.inv(self._lattice)
        self._positions = np.dot(self._lattice, self._supercell.get_points())

    def _get_modulation(self, point):
        modulation = self._add_modulations(point)
        modulation *= self._get_normalize_phase_factor(modulation)
        return modulation
        
    def _add_modulations(self, phase_set):
        modulation = np.zeros(self._vectors[0].shape, dtype=complex)
        for c, v in zip(phase_set, self._vectors):
            modulation += c * v
        return modulation
        
    def _get_cell_with_modulation(self, modulation):
        points = np.dot(self._lattice_inv,
                        self._positions + modulation.real)
        for p in points.T:
            p -= np.floor(p)
    
        return Cell(lattice=self._lattice,
                    points=points,
                    masses=self._supercell.get_masses(),
                    numbers=self._supercell.get_numbers())
    
    def _get_normalize_phase_factor(self, modulation):
        u = modulation.flatten()
        index_max_elem = np.argmax(abs(u))
        max_elem = u[index_max_elem]
        return max_elem / abs(max_elem)

    def _get_normalize_amplitude(self, modulation):
        u = modulation.flatten().real
        index_max_elem = np.argmax(abs(u))
        max_elem = u[index_max_elem]
        return self._max_displacement / abs(max_elem)

    def _get_phases(self):
        n = len(self._vectors)
        phases = np.exp(np.arange(self._ndiv) * 2j * np.pi / self._ndiv)

        if n == 1:
            phase_sets = [[1.]]
        elif n == 2:
            phase_sets = []
            for x in phases:
                phase_sets.append([1, x])
        elif n == 3:
            phase_sets = []
            for x in phases:
                for y in phases:
                    phase_sets.append([1, x, y])
        elif n == 4:
            phase_sets = []
            for x in phases:
                for y in phases:
                    for z in phases:
                        phase_sets.append([1, x, y, z])
        elif n == 5:
            phase_sets = []
            for x_1 in phases:
                for x_2 in phases:
                    for x_3 in phases:
                        for x_4 in phases:
                            phase_sets.append([1, x_1, x_2, x_3, x_4])
        elif n == 6:
            phase_sets = []
            for x_1 in phases:
                for x_2 in phases:
                    for x_3 in phases:
                        for x_4 in phases:
                            for x_5 in phases:
                                phase_sets.append([1, x_1, x_2, x_3, x_4, x_5])


        return phase_sets

    def _get_phase_shifts_at_lattice_points(self):
        total_phase_shifts = [1+0j]
        dim = np.array(self._modulation_dimension) * 3
        for i in range(self._modulation_dimension[0]):
            for j in range(self._modulation_dimension[1]):
                for k in range(self._modulation_dimension[2]):
                    phase = np.exp(2j * np.pi * 
                                   np.dot([i, j, k], 1. / dim))
                    is_found = True
                    for p in total_phase_shifts:
                        if abs(p - phase) < 1e-10:
                            is_found = False
                            break
                    if is_found:
                        total_phase_shifts.append(phase)

        return total_phase_shifts
        

class PhononModulationOld:
    def __init__(self,
                 phonon,
                 qpoint,
                 band_indices,
                 modulation_dimension,
                 ndiv=180,
                 symmetry_tolerance=0.05,
                 max_displacement=0.1,
                 store_all=False):
        self._phonon = phonon
        self._qpoint = qpoint
        self._band_indices = band_indices
        self._modulation_dimension = modulation_dimension
        self._ndiv = ndiv
        self._symmetry_tolerance = symmetry_tolerance
        self._max_displacement = max_displacement
        self._store_all = store_all
        self._vectors = []
        self._arguments = []
        self._points_on_sphere = []
        self._all_cells = []

        self._run()
        
    def get_modulation_cells(self):
        return [self._get_cell_with_modulation(m)
                for m in self.get_modulations()]

    def get_modulations(self):
        return [self._get_modulation(p) for p in self._points_on_sphere]

    def get_all_cells(self):
        return self._all_cells

    def get_points_on_sphere(self):
        return self._points_on_sphere

    def get_vectors(self):
        return self._vectors

    def get_arguments(self):
        return self._arguments

    def get_supercell(self):
        return self._supercell

    def _run(self):
        self._set_best_arguments_of_vectors_and_supercell()
        max_num_op = 0
        best_cells = []
        points_on_sphere = []
        for i, point in enumerate(self._get_all_points_on_sphere()):
            modcell = self._get_cell_with_modulation(
                self._get_modulation(point))
            symmetry = get_symmetry_dataset(
                modcell, tolerance=self._symmetry_tolerance)

            if self._store_all:
                self._all_cells.append(modcell)

            refined_cell = get_crystallographic_cell(
                modcell,
                tolerance=self._symmetry_tolerance)
            num_op = len(symmetry['rotations'])
            if num_op > max_num_op:
                max_num_op = num_op
                best_cells = [refined_cell]
                points_on_sphere = [point.copy()]
            if num_op == max_num_op:
                is_found = True
                for bc in best_cells:
                    if xtal_compare(bc,
                                    refined_cell,
                                    tolerance=self._symmetry_tolerance,
                                    angle_tolerance=1.0):
                        is_found = False
                        break
                if is_found:
                    best_cells.append(refined_cell)
                    points_on_sphere.append(point.copy())

        self._points_on_sphere = points_on_sphere

    def _set_best_arguments_of_vectors_and_supercell(self):
        phonon_modes = [[self._qpoint, i, 1, 0]
                        for i in self._band_indices]
        self._phonon.set_modulations(self._modulation_dimension, phonon_modes)
        modulations, supercell = self._phonon.get_delta_modulations()
            
        self._supercell = atoms2cell(supercell)

        vectors = []
        arguments = []
        for i, deltas in enumerate(modulations):
            argument, indices = self._best_argument(deltas.T)
            arguments.append([argument, indices])
            vectors.append((deltas.T * np.exp(-1j * argument)).real)

        self._vectors = vectors
        self._arguments = arguments
    
    def _get_modulation(self, point):
        modulation = np.zeros(self._vectors[0].shape, dtype=float)
        for c, v in zip(point, self._vectors):
            modulation += c * v
        return modulation
        
    def _get_all_points_on_sphere(self):
        ndiv = self._ndiv
        n = len(self._vectors)
        phi = np.linspace(0, 2 * np.pi, ndiv * 2)
        theta = np.linspace(0, np.pi, ndiv)
    
        if n == 1:
            x_i = np.array([[1.0]])
    
        elif n == 2:
            x_i = np.transpose([np.cos(phi), np.sin(phi)])
    
        elif n == 3:
            x = []
            for a in np.sin(theta):
                for b in np.cos(phi):
                    x.append(a * b)
            y = []
            for a in np.sin(theta):
                for b in np.sin(phi):
                    y.append(a * b)
            z = []
            for a in np.cos(theta):
                for b in np.cos(phi):
                    z.append(a * 1)
            x_i = np.transpose([x, y, z])
    
        elif n == 4:
            x1 = []
            for a in np.cos(theta): # C
                for b in range(ndiv):
                    for c in range(ndiv * 2):
                        x1.append(a * 1 * 1)
            x2 = []
            for a in np.sin(theta): # SC
                for b in np.cos(theta):
                    for c in range(ndiv * 2):
                        x2.append(a * b * 1)
            x3 = []
            for a in np.sin(theta): # SSC
                for b in np.sin(theta):
                    for c in np.cos(phi):
                        x3.append(a * b * c)
            x4 = []
            for a in np.sin(theta): # SSS
                for b in np.sin(theta):
                    for c in np.sin(phi):
                        x4.append(a * b * c)
            x_i = np.transpose([x1, x2, x3, x4])
    
        elif n == 6:
            x1 = []
            for a in np.cos(theta): # C
                for b in range(ndiv):
                    for c in range(ndiv):
                        for d in range(ndiv):
                            for e in range(ndiv * 2):
                                x1.append(a * 1 * 1 * 1 * 1)
            x2 = []
            for a in np.sin(theta): # SC
                for b in np.cos(theta):
                    for c in range(ndiv):
                        for d in range(ndiv):
                            for e in range(ndiv * 2):
                                x2.append(a * b * 1 * 1 * 1)
            x3 = []
            for a in np.sin(theta): # SSC
                for b in np.sin(theta):
                    for c in np.cos(theta):
                        for d in range(ndiv):
                            for e in range(ndiv * 2):
                                x3.append(a * b * c * 1 * 1)
            x4 = []
            for a in np.sin(theta): # SSSC
                for b in np.sin(theta):
                    for c in np.sin(theta):
                        for d in np.cos(theta):
                            for e in range(ndiv * 2):
                                x4.append(a * b * c * d * 1)
            x5 = []
            for a in np.sin(theta): # SSSSC
                for b in np.sin(theta):
                    for c in np.sin(theta):
                        for d in np.sin(theta):
                            for e in np.cos(phi):
                                x5.append(a * b * c * d * e)
            x6 = []
            for a in np.sin(theta): # SSSSS
                for b in np.sin(theta):
                    for c in np.sin(theta):
                        for d in np.sin(theta):
                            for e in np.sin(phi):
                                x6.append(a * b * c * d * e)
    
    
            x_i = np.transpose([x1, x2, x3, x4, x5, x6])
    
        else:
            x_i = None
                        
        return x_i
    
    def _get_cell_with_modulation(self, modulation):
        supercell = self._supercell
        lattice = supercell.get_lattice()
        max_mod = 0.0
        for mod in modulation.T:
            if max_mod < np.linalg.norm(mod):
                max_mod = np.linalg.norm(mod)

        points = np.dot(np.linalg.inv(lattice),
                        np.dot(lattice, supercell.get_points()) +
                        modulation / max_mod * self._max_displacement)
        for p in points.T:
            p -= np.floor(p)
    
        return Cell(lattice=lattice,
                    points=points,
                    masses=supercell.get_masses(),
                    numbers=supercell.get_numbers())
    
    def _best_argument(self, deltas):
        num_atom = deltas.shape[1] / np.prod(self._modulation_dimension)
        argument = 0
        max_val = 0
        best_indices = None
        for i in range(num_atom):
            for j in range(3):
                if abs(deltas[j, i]) > max_val:
                    max_val = abs(deltas[j, i])
                    argument = np.angle(deltas[j, i])
                    best_indices = (j, i)
    
        return argument, best_indices
