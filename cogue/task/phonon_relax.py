"""Relax crystal structure constrained by symmetry

1. Structure optimization

2. Calculate force constants (second order coefficients of crystal
   potential with respect to atomic displacements) of some specific
   size of supercell.

3. Find negative eigenvalues of the force constant matrix, then
   displace atomic coordinates along the corresponding
   eigenvectors. This deformation breaks the crystal symmetry. In this
   task, phonon eigenvector, which is eigenvector of dynamical matrix,
   is employed as the eigenvector.

4. Structure optimization with tight residual force convergence criterion.

5. Search the crystal symmetry after the precise relaxation. Higher
   symmetry than that of distorted structure is expected if possible.

"""

import os
import numpy as np
from cogue.task import TaskElement
from cogue.crystal.cell import Cell
from cogue.crystal.converter import \
    atoms2cell, write_cif_P1, write_v_sim, get_lattice_parameters
from cogue.calculator.vasp.vasp_io import write_poscar
from cogue.crystal.symmetry import \
    get_symmetry_dataset, get_crystallographic_cell
from cogue.crystal.converter import get_primitive
from cogue.crystal.xtalcomp import compare as xtal_compare

CUTOFF_ZERO = 1e-10
DEGENERACY_TOLERANCE = 1e-3

class PhononRelaxBase(TaskElement):
    def __init__(self,
                 directory=None,
                 name=None,
                 ancestral_cells=[],
                 distance=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 symmetry_tolerance=None,
                 restrict_offspring=False,
                 max_offspring=None,
                 cutoff_eigenvalue=None,
                 max_displacement=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._ancestral_cells = ancestral_cells
        self._task_type = "phonon_relax"
        self._distance = distance
        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._max_iteration = max_iteration
        self._min_iteration = min_iteration
        self._symmetry_tolerance = symmetry_tolerance
        self._restrict_offspring = restrict_offspring
        self._max_offspring = max_offspring
        self._cutoff_eigenvalue = cutoff_eigenvalue
        if max_displacement:
            self._max_displacement = max_displacement
        else:
            self._max_displacement = symmetry_tolerance
        self._traverse = traverse

        self._phr_tasks = []
        self._tasks = []
        self._stage = 0

        self._energy = None

    def set_status(self):
        done = True
        terminate = False
        reversion = False

        if self._stage == 0 and self._tasks[0].get_status() == "reversion":
            reversion = True

        if not reversion:
            for task in self._tasks:
                done &= task.done()
                if task.get_status() == "terminate":
                    terminate = True

        if done:
            if terminate:
                self._status = "terminate"
            elif reversion:
                self._status = "reversion"
            else:
                self._status = "next"

        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        self._overwrite_settings()

        self._status = "stage 0"
        self._stage = 0
        self._tasks = []
        task = self._get_phonon_relax_element_task(self._cell)
        self._phr_tasks = [task]
        self._tasks = [task]
        space_group = get_symmetry_dataset(self._cell,
                                           tolerance=self._symmetry_tolerance)
        self._comment = space_group['international']

    def end(self):
        self._write_yaml()

    def done(self):
        return (self._status == "terminate" or
                self._status == "reversion" or
                self._status == "next" or
                self._status == "done")

    def next(self):
        if self._stage == 0:
            task = self._phr_tasks[0]
            if self._status == "next":
                imag_modes = task.get_imaginary_modes()
                self._comment += " --> %s" % task.get_space_group_type()
                self._energy = task.get_energy()
                if self._energy:
                    self._comment += "\\n%f" % self._energy
                if imag_modes: # Next phonon-relaxes
                    self._stage = 1
                    self._status = "offspring"
                    self._set_offsprings(imag_modes)
                    return self._tasks
                else: # No imaginary mode
                    self._status = "done"
                    self._tasks = []
                    raise StopIteration
            elif self._status == "terminate":
                raise StopIteration
            elif self._status == "reversion":
                self._comment += " --> %s" % task.get_space_group_type()
                self._tasks = []
                raise StopIteration
            else:
                print "It is something wrong in phonon_relax.py."
                raise StopIteration
        else:
            if self._status == "next":
                self._status = "done"
            raise StopIteration

    def _set_offsprings(self, imag_modes):
        self._tasks = []
        # If self._restrict_offspring == True and
        # number of imaginary modes at Gamma > 1,
        # offspring is restricted only at Gamma.
        # If self._restrict_offspring is a positive integer,
        # the restriction is set as number of imaginary modes
        # at Gamma > self._restrict_offspring.
        num_gamma = [x[3] for x in imag_modes].count(0)
        if (self._restrict_offspring and 
            num_gamma > self._restrict_offspring * 1):
            mod_modes = []
            for i, x in enumerate(imag_modes):
                if x[3] == 0:
                    mod_modes.append((x[0], x[2], i)) # cell & Im(freq)
        else:
            mod_modes = [(x[0], x[2], i) for i, x in enumerate(imag_modes)]
        mod_modes.sort(key=lambda mod_modes: -mod_modes[1])
        if self._max_offspring:
            mod_modes = mod_modes[:min(self._max_offspring, len(mod_modes))]
        for (cell, freq, index) in mod_modes:
            self._tasks.append(self._get_phonon_relax_task(
                    cell,
                    self._ancestral_cells,
                    "phonon_relax-%d" % (index + 1)))
        self._phr_tasks += self._tasks
        
    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
        w.write("pressure_target: %f\n" % self._pressure_target)
        w.write("stress_tolerance: %f\n" % self._stress_tolerance)
        w.write("force_tolerance: %f\n" % self._force_tolerance)
        w.write("max_increase: %f\n" % self._max_increase)
        w.write("max_iteration: %d\n" % self._max_iteration)
        w.write("min_iteration: %d\n" % self._min_iteration)
        w.write("symmetry_tolerance: %d\n" % self._symmetry_tolerance)
        if self._restrict_offspring:
            w.write("restrict_offspring: %s\n" % self._restrict_offspring)
        if self._max_offspring:
            w.write("max_offspring: %s\n" % self._max_offspring)
        if self._energy:
            w.write("electric_total_energy: %20.10f\n" % self._energy)
        w.write("status: %s\n" % self._status)
        w.write("tasks:\n")

        for task in self._phr_tasks:
            if task.get_status():
                w.write("- name:   %s\n" % task.get_name())
                w.write("  status: %s\n" % task.get_status())
        w.close()

class PhononRelaxElementBase(TaskElement):
    """PhononRelaxElementBase class

    Relax crystal structure using phonon eigenvector with imaginary
    frequency

    """

    def __init__(self,
                 directory=None,
                 name=None,
                 ancestral_cells=[],
                 distance=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 symmetry_tolerance=None,
                 cutoff_eigenvalue=None,
                 max_displacement=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._ancestral_cells = ancestral_cells
        self._task_type = "phonon_relax_element"
        self._distance = distance
        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._max_iteration = max_iteration
        self._min_iteration = min_iteration
        self._symmetry_tolerance = symmetry_tolerance
        self._cutoff_eigenvalue = cutoff_eigenvalue
        self._max_displacement = max_displacement
        self._traverse = traverse

        self._tasks = []
        self._phre = []
        self._stage = 0
        self._supercell_dimention = None
        self._space_group_type = None
        self._energy = None
        # self._imaginary_modes:
        # (cell with modulation,
        #  q-point,
        #  Im(frequency),
        #  Index of q-point in the mesh sampling,
        #  Band index,
        #  Number of degeneracy)
        self._imaginary_modes = []

    def get_imaginary_modes(self):
        return self._imaginary_modes

    def get_space_group_type(self):
        return self._space_group_type

    def get_energy(self):
        return self._energy

    def set_status(self):
        all_done = True
        for task in self._tasks:
            status = task.get_status()
            all_done &= (status == "done" or status == "max_iteration")
            if status == "terminate":
                self._status = "terminate"
                break
            
        if all_done:
            self._status = "next"

        if self._space_group_type:
            self._comment = self._space_group_type
            
        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        self._overwrite_settings()

        self._status = "stage 0"
        self._stage = 0
        self._set_stage0()

    def end(self):
        self._comment = self._space_group_type
        self._write_yaml()

    def done(self):
        return ("terminate" in self._status or
                "reversion" in self._status or 
                "done" in self._status or
                "next" in self._status)

    def next(self):
        if self._stage == 0:
            if self._status == "terminate":
                raise StopIteration
            else:
                cell = self._tasks[0].get_cell()
                for ancest in self._ancestral_cells:
                    if xtal_compare(ancest,
                                    cell,
                                    tolerance=self._symmetry_tolerance):
                        symmetry = get_symmetry_dataset(
                            cell,
                            tolerance=self._symmetry_tolerance)
                        self._space_group_type = symmetry['international']
                        self._status = "reversion"
                        raise StopIteration
                self._ancestral_cells.append(cell)
                self._set_stage1(cell)
                self._stage = 1
                self._status = "stage 1"
                return self._tasks
        else:
            if self._status == "next":
                self._analyze_phonon()
                self._status = "done"

            raise StopIteration

    def _set_stage0(self):
        self._tasks = []
        task = self._get_equilibrium_task(
            cell=self._cell,
            impose_symmetry=True,
            symmetry_tolerance=self._symmetry_tolerance)
        symmetry = get_symmetry_dataset(self._cell,
                                        tolerance=self._symmetry_tolerance)
        self._space_group_type = symmetry['international']
        self._phre_tasks = [task]
        self._tasks = [task]

    def _set_stage1(self, cell):
        prim_cell = get_primitive(cell, tolerance=self._symmetry_tolerance)
        sym_dataset = get_symmetry_dataset(prim_cell)
        self._space_group_type = sym_dataset['international']
        spg_number = sym_dataset['number']
        if (spg_number >= 143 and
            spg_number <= 194 and
            not self._space_group_type[0] == 'R'): # Hexagonal lattice
            self._supercell_dimensions = [[3, 3, 2], [2, 2, 2]]
        else: # Other cases
            self._supercell_dimensions = [[2, 2, 2]]

        # Long cell axis is not multiplied.
        for dimension in self._supercell_dimensions:
            for i, length in enumerate(
                get_lattice_parameters(prim_cell.get_lattice())):
                if length * dimension[i] > 20:
                    dimension[i] = 1

        self._tasks = []
        for i, dimension in enumerate(self._supercell_dimensions):
            task = self._get_phonon_task(prim_cell,
                                         np.diag(dimension),
                                         "phonon-%d" % (i + 1))
            self._phre_tasks.append(task)
            self._tasks.append(task)

    def _analyze_phonon(self):
        for dimension, task in zip(self._supercell_dimensions,
                                   self._tasks):
            self._energy = task.get_energy()
            phonon = task.get_phonon()
            phonon.set_mesh(dimension, is_gamma_center=True)
            qpoints, weigths, frequencies, eigvecs = phonon.get_mesh()
            eigenvalues = frequencies ** 2 * np.sign(frequencies)
            if (eigenvalues < self._cutoff_eigenvalue).any():
                imaginary_modes = []
                qpoints_done = [imag_mode[1]
                                for imag_mode in self._imaginary_modes]
                self._imaginary_modes += get_unstable_modulations(
                    phonon,
                    dimension,
                    max_displacement=self._max_displacement,
                    cutoff_eigenvalue=self._cutoff_eigenvalue,
                    ndiv=180,
                    excluded_qpoints=qpoints_done)
                
    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
        w.write("pressure_target: %f\n" % self._pressure_target)
        w.write("stress_tolerance: %f\n" % self._stress_tolerance)
        w.write("force_tolerance: %f\n" % self._force_tolerance)
        w.write("max_increase: %f\n" % self._max_increase)
        w.write("max_iteration: %d\n" % self._max_iteration)
        w.write("min_iteration: %d\n" % self._min_iteration)
        w.write("symmetry_tolerance: %f\n" % self._symmetry_tolerance)
        w.write("stage: %d\n" % self._stage)
        w.write("status: %s\n" % self._status)
        if self._energy:
            w.write("electric_total_energy: %20.10f\n" % self._energy)
        if self._imaginary_modes:
            w.write("imaginary_modes:\n")
            for imag_mode in self._imaginary_modes:
                spg = get_symmetry_dataset(imag_mode[0],
                                           tolerance=self._symmetry_tolerance)
                q = imag_mode[1]
                freq = imag_mode[2]
                q_index = imag_mode[3] + 1
                band_index = imag_mode[4] + 1
                degeneracy = imag_mode[6]
                dimension = tuple(imag_mode[7])
                w.write("- supercell_dimension: [ %d, %d, %d ]\n" % dimension)
                w.write("  qpoint: [ %6.4f, %6.4f, %6.4f ] # %d\n" %
                        (q[0], q[1], q[2], q_index))
                w.write("  frequency: %10.5f # %d\n" % (-freq, band_index))
                w.write("  degeneracy: %d\n" % degeneracy)
                w.write("  space_group_type: %s\n" % spg['international'])
                w.write("  space_group_number: %s\n" % spg['number'])
        w.write("tasks:\n")
        for task in self._phre_tasks:
            if task.get_status():
                w.write("- name:   %s\n" % task.get_name())
                w.write("  status: %s\n" % task.get_status())
        w.close()


def get_unstable_modulations(phonon,
                             supercell_dimension,
                             degeneracy_tolerance=DEGENERACY_TOLERANCE,
                             max_displacement=None,
                             cutoff_eigenvalue=None,
                             ndiv=180,
                             excluded_qpoints=None):
    qpoints, weigths, frequencies, eigvecs = phonon.get_mesh()
    eigenvalues = frequencies ** 2 * np.sign(frequencies)
    imag_modes = []

    for i, (q, eigs_at_q) in enumerate(zip(qpoints, eigenvalues)):
        qpt_exists = False
        for qpt in excluded_qpoints:
            if (abs(q - qpt) < 1e-10).all():
                qpt_exists = True
                break
        if qpt_exists:
            continue
        
        indices_imaginary = np.where(eigs_at_q < cutoff_eigenvalue)[0]
        degeneracy_sets = get_degeneracy_sets(eigs_at_q,
                                              indices_imaginary,
                                              degeneracy_tolerance)
        if degeneracy_sets:
            phonon.write_animation(q, filename="anime-d%d%d%d-q%d.ascii" %
                                   (tuple(supercell_dimension) + (i + 1,)))
        
        for deg_set in degeneracy_sets:
            j = deg_set[0]
            eig = eigs_at_q[j]
            modulation_dimension = []
            for a, multi in zip(q, supercell_dimension):
                if abs(a) < CUTOFF_ZERO:
                    modulation_dimension.append(1)
                else:
                    modulation_dimension.append(multi)

            phononMod = PhononModulation(
                phonon,
                q,
                deg_set,
                modulation_dimension,
                ndiv=ndiv,
                symmetry_tolerance=max_displacement * 0.9,
                max_displacement=max_displacement,
                store_all=False)
            # write_all_symmetries(phononMod.get_all_symmetries(), i , j)
            modulation_cells = phononMod.get_modulation_cells()
            supercell = phononMod.get_supercell()
            write_poscar(supercell,
                         "SPOSCAR-d%d%d%d" % tuple(modulation_dimension))
            write_cif_P1(supercell,
                         "supercell-d%d%d%d.cif" % tuple(modulation_dimension))
            write_v_sim(supercell,
                        "supercell-d%d%d%d.ascii" % tuple(modulation_dimension))
            for k, modcell in enumerate(modulation_cells):
                write_poscar(modcell,
                             "POSCAR-d%d%d%d-q%db%d-%d" %
                             (tuple(modulation_dimension) + (i + 1, j + 1, k + 1)))
                write_cif_P1(modcell,
                             "mod-d%d%d%d-q%db%d-%d.cif" % 
                             (tuple(modulation_dimension) + (i + 1, j + 1, k + 1)))
                write_v_sim(modcell,
                            "mod-d%d%d%d-q%db%d-%d.ascii" % 
                            (tuple(modulation_dimension) + (i + 1, j + 1, k + 1)))
                imag_modes.append(
                    (modcell, q, np.sqrt(-eig), i, j, k,
                     len(deg_set), modulation_dimension))
            
    return imag_modes

def write_all_symmetries(all_symmetries, i, j):
    """ Write all symmetry information

    To obtain all_symmetries, PhononModulation has to be invoked with
    store_all=True.
    
    """
    w = open("symmetries-q%db%d.dat" % (i + 1, j + 1), 'w')
    for sym in all_symmetries:
        w.write("%s %d\n" % (sym['international'],
                             len(sym['rotations'])))
    w.close()

def get_degeneracy_sets(eigs, indices, degeneracy_tolerance):
    degeneracy_sets = []
    for j in indices:
        deg_set = []
        for k in indices:
            if abs(eigs[j] - eigs[k]) < degeneracy_tolerance:
                deg_set.append(k)
        if j == deg_set[0]:
            degeneracy_sets.append(deg_set)

    return degeneracy_sets

class PhononModulation:
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
                                    tolerance=self._symmetry_tolerance):
                        is_found = False
                        break
                if is_found:
                    best_cells.append(refined_cell)
                    points_on_sphere.append(point.copy())

        self._points_on_sphere = points_on_sphere

    def _set_best_arguments_of_vectors_and_supercell(self):
        modulations, supercell = self._phonon.get_delta_modulations(
            self._qpoint, self._modulation_dimension)
        self._supercell = atoms2cell(supercell)

        vectors = []
        arguments = []
        for i, deltas in enumerate(modulations):
            if i in self._band_indices:
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
    



if __name__ == '__main__':
    import sys
    from phonopy import Phonopy
    from phonopy.structure.atoms import Atoms
    from phonopy.structure.symmetry import Symmetry
    from phonopy.interface.vasp import read_vasp, write_vasp
    from phonopy.hphonopy.file_IO import parse_FORCE_CONSTANTS
    from optparse import OptionParser
    
    def get_parameters():
        parser = OptionParser()
        parser.set_defaults(supercell_dimension=None,
                            qpoint=None,
                            band_indices=None,
                            ndiv=None,
                            tolerance=None)
        parser.add_option("--dim", dest="supercell_dimension",
                          action="store", type="string")
        parser.add_option("--moddim", dest="modulation_dimension",
                          action="store", type="string")
        parser.add_option("--ndiv", dest="ndiv",
                          action="store", type="int")
        parser.add_option("--band", dest="band_indices",
                          action="store", type="string")
        parser.add_option("-q", dest="qpoint",
                          action="store", type="string")
        parser.add_option("--tolerance", dest="tolerance", type="float",
                          help="Symmetry tolerance to search")
        (options, args) = parser.parse_args()
        
        if not options.supercell_dimension:
            print "Option --dim has to be set."
            sys.exit(1)
        
        if not options.qpoint:
            print "Option -q has to be set."
            sys.exit(1)
    
        if not options.band_indices:
            print "Option --band has to be set."
            sys.exit(1)
    
        qpoint = [float(x) for x in options.qpoint.split()]
        if len(qpoint) == 3:
            qpoint = np.array(qpoint)
        else:
            print "Illeagal q-point"
            sys.exit(1)
        
        supercell_dimension = \
            [int(x) for x in options.supercell_dimension.split()]
        if not len(supercell_dimension) == 3:
            print "Illeagal supercell dimension"
            sys.exit(1)
    
        modulation_dimension = \
            [int(x) for x in options.modulation_dimension.split()]
        if not len(modulation_dimension) == 3:
            print "Illeagal modulatoin dimension"
            sys.exit(1)
    
        band_indices = [int(x) - 1 for x in options.band_indices.split()]
        if len(band_indices) > 4:
            print "Number of band indices is too large (%d)" % len(band_indices)
            sys.exit(1)
    
        if options.ndiv:
            ndiv = options.ndiv
        else:
            ndiv = 360
    
        if options.tolerance:
            tolerance = options.tolerance
        else:
            tolerance = 1e-5
    
        cell = read_vasp(args[0])
        fc = parse_FORCE_CONSTANTS(args[1])
    
        return (qpoint,
                supercell_dimension,
                modulation_dimension,
                cell,
                fc,
                ndiv,
                band_indices,
                tolerance)

    def get_phonon(cell, supercell_dimension, fc):
        phonon = Phonopy(cell,
                         np.diag(supercell_dimension),
                         is_auto_displacements=False)
        phonon.set_post_process(force_constants=fc)
    
        return phonon


    (qpoint,
     supercell_dimension,
     modulation_dimension,
     cell,
     fc,
     ndiv,
     band_indices,
     tolerance) = get_parameters()

    phonon = get_phonon(cell, supercell_dimension, fc)
    phononMod = PhononModulation(phonon,
                                 qpoint,
                                 band_indices,
                                 modulation_dimension,
                                 ndiv=ndiv,
                                 symmetry_tolerance=tolerance,
                                 max_displacement=tolerance * 1.1)

    best_cells = phononMod.get_modulation_cells()
    for cell in best_cells:
        sym = get_symmetry_dataset(cell, tolerance)
        print sym['international'], len(sym['rotations'])
