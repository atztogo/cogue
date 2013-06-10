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
from cogue.crystal.converter import \
    write_cif_P1, write_v_sim, get_lattice_parameters
from cogue.interface.vasp_io import write_poscar
from cogue.crystal.converter import get_primitive
from cogue.crystal.xtalcomp import compare as xtal_compare
from cogue.phonon.modulation import PhononModulation
from cogue.crystal.symmetry import \
     get_symmetry_dataset, get_crystallographic_cell

CUTOFF_ZERO = 1e-10
DEGENERACY_TOLERANCE = 1e-3
MAX_DISPLACEMENT_RATIO = 1.1

class PhononRelaxBase(TaskElement):
    def __init__(self,
                 directory=None,
                 name=None,
                 ancestral_cells={},
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
                 num_sampling_points=None,
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
            self._max_displacement = symmetry_tolerance * MAX_DISPLACEMENT_RATIO
        self._num_sampling_points = num_sampling_points
        self._traverse = traverse

        self._phr_tasks = []
        self._tasks = []
        self._stage = 0

        self._energy = None

    def set_status(self):
        done = True
        terminate = False
        confluence = None

        if self._stage == 0 and "confluence" in self._tasks[0].get_status():
            confluence = self._tasks[0].get_status()

        if not confluence:
            for task in self._tasks:
                done &= task.done()
                if task.get_status() == "terminate":
                    terminate = True

        if done:
            if terminate:
                self._status = "terminate"
            elif confluence:
                self._status = confluence
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
        self._comment = space_group['international_standard']

    def end(self):
        self._write_yaml()

    def done(self):
        return ("terminate" in self._status or
                "confluence" in self._status or
                "next" in self._status or
                "done" in self._status)

    def next(self):
        if self._stage == 0:
            task = self._phr_tasks[0]
            if self._status == "next":
                imag_modes = task.get_imaginary_modes()
                self._comment += " --> %s" % task.get_space_group_type()
                self._energy = task.get_energy()
                if self._energy:
                    num_atom = len(task.get_cell().get_symbols())
                    self._comment += "\\n%f/%d" % (self._energy, num_atom)
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
            elif "confluence" in self._status:
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
        w.write("max_displacement: %f\n" % self._max_displacement)
        w.write("cutoff_eigenvalue: %f\n" % self._cutoff_eigenvalue)
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
                 ancestral_cells={},
                 tid_parent=None,
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
                 num_sampling_points=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._ancestral_cells = ancestral_cells
        self._tid_parent = tid_parent
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
        self._num_sampling_points = num_sampling_points
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

    def get_cell(self):
        if self._stage == 0:
            return self._cell
        else:
            return self._phre_tasks[1].get_cell()

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
                "confluence" in self._status or 
                "done" in self._status or
                "next" in self._status)

    def next(self):
        if self._stage == 0:
            if self._status == "terminate":
                raise StopIteration

            cell = self._tasks[0].get_cell()
            tid = self._find_equivalent_crystal_structure(cell)
            if tid > 0: # Equivalent structure found
                self._status = "confluence with [%d]" % tid
                symmetry = get_symmetry_dataset(
                    cell,
                    tolerance=self._symmetry_tolerance)
                self._space_group_type = symmetry['international_standard']
                raise StopIteration
            elif (self._traverse ==  "restart" and 
                  not os.path.exists("phonon-1")):
                # This condition means the structure optimization terminated 
                # by that equivalent crystal strucutre was found. However
                # in restart mode, the order to parse directory tree can
                # be different from that in run time. So inequivalent
                # crystal structure can be found different point.
                # This condition indicates there should be an equivalent
                # crystal structure somewhere else. In this case, this
                # 'next' does nothing (restarting stage0) and waits for
                # until it will be found.
                self.begin()
                return self._tasks
            else: # No equivalent structure found, move to phonon calculation
                self._ancestral_cells[self._tid_parent] = cell
                self._set_stage1(cell)
                self._stage = 1
                self._status = "stage 1"
                return self._tasks
        else:
            if self._status == "next":
                self._analyze_phonon()
                self._status = "done"

            raise StopIteration

    def _find_equivalent_crystal_structure(self, cell):
        for tid in self._ancestral_cells:
            if xtal_compare(self._ancestral_cells[tid],
                            cell,
                            tolerance=self._symmetry_tolerance,
                            angle_tolerance=1.0):
                return tid
        return 0

    def _set_stage0(self):
        self._tasks = []
        task = self._get_equilibrium_task(
            cell=self._cell,
            impose_symmetry=True,
            symmetry_tolerance=self._symmetry_tolerance)
        symmetry = get_symmetry_dataset(self._cell,
                                        tolerance=self._symmetry_tolerance)
        self._space_group_type = symmetry['international_standard']
        self._phre_tasks = [task]
        self._tasks = [task]

    def _set_stage1(self, cell):
        prim_cell = get_primitive(cell, tolerance=self._symmetry_tolerance)
        sym_dataset = get_symmetry_dataset(prim_cell)
        self._space_group_type = sym_dataset['international_standard']
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
                    symmetry_tolerance=self._symmetry_tolerance,
                    max_displacement=self._max_displacement,
                    cutoff_eigenvalue=self._cutoff_eigenvalue,
                    ndiv=self._num_sampling_points,
                    excluded_qpoints=qpoints_done)

        sym_dataset = get_symmetry_dataset(
            self._tasks[0].get_cell())
        self._space_group_type = sym_dataset['international_standard']
                
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
        w.write("max_displacement: %f\n" % self._max_displacement)
        w.write("cutoff_eigenvalue: %f\n" % self._cutoff_eigenvalue)
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
                w.write("  band: %d\n" % band_index)
                w.write("  frequency: %10.5f\n" % (-freq))
                w.write("  degeneracy: %d\n" % degeneracy)
                w.write("  space_group_type: %s\n" % spg['international_standard'])
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
                             symmetry_tolerance=0.1,
                             max_displacement=0.2,
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
                symmetry_tolerance=symmetry_tolerance,
                max_displacement=max_displacement)
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

    



