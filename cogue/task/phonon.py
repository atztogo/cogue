import sys
import numpy as np
from cogue.task import TaskElement
from cogue.task.structure_optimization import StructureOptimizationYaml
from cogue.interface.vasp_io import write_poscar, write_poscar_yaml
from cogue.crystal.cell import sort_cell_by_symbols
from cogue.crystal.converter import cell2atoms
from cogue.crystal.supercell import estimate_supercell_matrix
from cogue.crystal.symmetry import get_crystallographic_cell

try:
    from phonopy import Phonopy
    from phonopy.interface import PhonopyYaml
    from phonopy.file_IO import write_FORCE_SETS
except ImportError:
    print("You need to install phonopy.")
    sys.exit(1)


class PhononYaml(StructureOptimizationYaml):
    def _get_phonon_yaml_lines(self, cell):
        lines = self._get_structopt_yaml_lines()
        lines.append("supercell_matrix:")
        if self._supercell_matrix is None:
            lines.append("- automatic")
        else:
            for row in self._supercell_matrix:
                lines.append("- [ %3d, %3d, %3d ]" % tuple(row))
        if self._primitive_matrix is not None:
            lines.append("primitive_matrix:")
            for row in self._primitive_matrix:
                lines.append("- [ %6.3f, %6.3f, %6.3f ]" % tuple(row))
        lines.append("distance: %f" % self._distance)
        if cell:
            lines += cell.get_yaml_lines()

        return lines


class PhononBase(TaskElement, PhononYaml):
    """PhononBase class

    This is an interface to phonopy.

    """

    def __init__(self,
                 directory=None,
                 name=None,
                 supercell_matrix=None,
                 primitive_matrix=None,
                 nac=False,
                 distance=None,
                 displace_plusminus='auto',
                 displace_diagonal=False,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 is_cell_relaxed=False,
                 max_num_atoms=None,
                 impose_symmetry=False,
                 stop_condition=None,
                 symmetry_tolerance=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "phonon"
        self._supercell_matrix = supercell_matrix
        if self._supercell_matrix is None:
            self._primitive_matrix = np.eye(3, dtype='double')
        else:
            self._primitive_matrix = primitive_matrix
        self._nac = nac
        self._distance = distance
        self._displace_plusminus = displace_plusminus
        self._displace_diagonal = displace_diagonal
        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._max_iteration = max_iteration
        self._min_iteration = min_iteration
        self._is_cell_relaxed = is_cell_relaxed
        self._max_num_atoms = max_num_atoms
        self._impose_symmetry = impose_symmetry
        self._stop_condition = stop_condition
        self._symmetry_tolerance = symmetry_tolerance
        self._traverse = traverse

        self._stage = 0
        self._tasks = []

        self._energy = None
        self._born = None
        self._epsilon = None

        self._space_group = None
        self._cell = None
        self._phonon = None  # Phonopy object
        self._all_tasks = None

        self._try_collect_forces = True

    def get_phonon(self):
        return self._phonon

    def get_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._all_tasks[0].get_cell()

    def get_energy(self):
        """Return energies at geometry optimization steps"""
        return self._energy

    def get_space_group(self):
        return self._space_group

    def set_status(self):
        if self._stage == 0:
            task = self._tasks[0]
            if task.done():
                self._space_group = task.get_space_group()
                status = task.get_status()
                if status == "done":
                    if not self._evaluate_stop_condition():
                        self._status = "next"
                else:
                    self._status = status
        else:
            done = True
            terminate = False
            for i, task in enumerate(self._tasks):
                done &= task.done()
                if (task.get_status() == "terminate" or
                    task.get_status() == "max_iteration"):
                    terminate = True
            if done:
                if terminate:
                    self._status = "terminate"
                else:
                    self._status = "next"

        self._write_yaml()

    def begin(self):
        if not self._job:
            print("set_job has to be executed.")
            raise RuntimeError

        if self._is_cell_relaxed:
            self._all_tasks = [None]
            self._set_stage1()
        else:
            self._set_stage0()

    def done(self):
        return (self._status == "terminate" or
                self._status == "done" or
                self._status == "max_iteration" or
                self._status == "low_symmetry" or
                self._status == "force_collection_failure" or
                self._status == "next")

    def __next__(self):
        return self.next()

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._energy = self._tasks[0].get_energy()
                num_atom = len(self._tasks[0].get_cell().get_symbols())
                self._comment = self._space_group['international']
                self._comment += "\\n%f/%d" % (self._energy, num_atom)
                self._set_stage1()
                return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._set_stage0()
                return self._tasks

        elif self._stage == 1:  # task 1..n: displaced supercells
            if self._status == "next":
                if self._collect_forces():
                    if self._nac:
                        self._set_stage2()
                        return self._tasks
                    else:
                        self._status = "done"
                else:
                    if self._try_collect_forces:
                        self._status = "displacements"
                        self._log += ("Collection of forces failed. "
                                      "Try once more.\n")
                        self._try_collect_forces = False
                        raise StopIteration
                    else:
                        self._status = "force_collection_failure"
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                disp_terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        disp_terminated.append(i)
                tasks = self._get_displacement_tasks()
                self._tasks = []
                for i in disp_terminated:
                    self._tasks.append(tasks[i])
                    self._all_tasks[i + 1] = tasks[i]
                self._status = "displacements"
                return self._tasks
        elif self._stage == 2:
            if self._status == "next":
                self._status = "done"
                self._set_born_and_epsilon()
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                self._all_tasks.pop()
                self._set_stage2()
        else:
            pass

        self._tasks = []
        self._write_yaml()
        raise StopIteration

    def _set_stage0(self):
        self._status = "equilibrium"
        task = self._get_equilibrium_task(
            impose_symmetry=self._impose_symmetry)
        self._all_tasks = [task]
        self._tasks = [task]

    def _set_stage1(self):
        self._stage = 1
        self._status = "displacements"
        self._set_phonon()
        self._tasks = self._get_displacement_tasks()
        self._all_tasks += self._tasks

    def _set_stage2(self):
        self._stage = 2
        self._status = "nac"
        if self._nac == "relax":
            nac_task = self._get_nac_task(is_cell_relaxed=False)
        else:
            nac_task = self._get_nac_task()
        self._tasks = [nac_task]
        self._all_tasks += self._tasks

    def _collect_forces(self):
        forces = []
        for task in self._tasks:
            forces.append(task.get_properties()['forces'][-1])
        if self._phonon.produce_force_constants(forces=forces):
            write_FORCE_SETS(self._phonon.get_displacement_dataset())
            return True
        else:
            # This can be due to delay of writting file to file system.
            return False

    def _set_born_and_epsilon(self):
        nac_task = self._tasks[0]
        born = nac_task.get_born_effective_charge()
        epsilon = nac_task.get_dielectric_constant()

        indep_atoms = self._phonon.get_symmetry().get_independent_atoms()
        supercell = self._phonon.get_supercell()
        s2u = supercell.get_supercell_to_unitcell_map()
        u2u = supercell.get_unitcell_to_unitcell_map()
        indep_atoms_u = [u2u[i] for i in s2u[indep_atoms]]

        if born is not None and epsilon is not None:
            self._born = born
            self._epsilon = epsilon
            header = "# epsilon and Z* of atoms "
            header += ' '.join(["%d" % (n + 1) for n in indep_atoms_u])
            lines = [header]
            lines.append(("%13.8f" * 9) % tuple(epsilon.flatten()))
            for z in born[indep_atoms_u]:
                lines.append(("%13.8f" * 9) % tuple(z.flatten()))
            with open("BORN", 'w') as w:
                w.write('\n'.join(lines))

    def _evaluate_stop_condition(self):
        if self._stop_condition:
            if "symmetry_operations" in self._stop_condition:
                num_ops = len(self._space_group['rotations'])
                if num_ops < self._stop_condition['symmetry_operations']:
                    self._status = "low_symmetry"
                    return True

        return False

    def _set_phonon(self):
        if self._supercell_matrix is None:
            cell = sort_cell_by_symbols(
                get_crystallographic_cell(self.get_cell(),
                                          tolerance=self._symmetry_tolerance))
            self._supercell_matrix = estimate_supercell_matrix(
                cell,
                max_num_atoms=self._max_num_atoms)
        else:
            cell = self.get_cell()

        phonopy_cell = cell2atoms(cell)
        self._phonon = Phonopy(phonopy_cell,
                               self._supercell_matrix,
                               primitive_matrix=self._primitive_matrix,
                               dynamical_matrix_decimals=14,
                               force_constants_decimals=14,
                               symprec=self._symmetry_tolerance)
        self._phonon.generate_displacements(
            distance=self._distance,
            is_plusminus=self._displace_plusminus,
            is_diagonal=self._displace_diagonal)

        write_poscar(cell, filename="POSCAR-unitcell")
        write_poscar_yaml(cell, filename="POSCAR-unitcell.yaml")
        phpy_yaml = PhonopyYaml(settings={'displacements': True})
        phpy_yaml.set_phonon_info(self._phonon)
        with open("phonopy_disp.yaml", 'w') as w:
            w.write(str(phpy_yaml))

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self.get_cell()
        lines += self._get_phonon_yaml_lines(cell)
        if self._all_tasks[0] is not None:
            if self._energy:
                lines.append("electric_total_energy: %20.10f" % self._energy)

        return lines
