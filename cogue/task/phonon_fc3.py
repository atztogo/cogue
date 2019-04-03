import sys
import numpy as np
from cogue.task import TaskElement
from cogue.interface.vasp_io import write_poscar
from cogue.crystal.converter import cell2atoms

try:
    from phonopy import Phonopy
except ImportError:
    print("You need to install phonopy.")
    sys.exit(1)
try:
    from phono3py.phonon3 import Phono3py
except ImportError:
    print("You need to install phono3py.")
    sys.exit(1)

from phonopy.file_IO import write_FORCE_SETS
from phonopy.file_IO import write_disp_yaml
from phono3py.file_IO import write_disp_fc3_yaml
from phono3py.file_IO import write_FORCES_FC3

class PhononFC3Base(TaskElement):
    """PhononFC3Base class

    This is an interface to anharmonic phonopy.

    """

    def __init__(self,
                 directory=None,
                 name=None,
                 supercell_matrix=None,
                 primitive_matrix=None,
                 with_perfect=True,
                 distance=None,
                 is_diagonal=True,
                 check_imaginary=True,
                 cutoff_frequency=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 is_cell_relaxed=False,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "anharmonic_phonon"
        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._with_perfect = with_perfect
        self._distance = distance
        self._is_diagonal = is_diagonal
        self._check_imaginary = check_imaginary
        self._cutoff_frequency = cutoff_frequency  # determine imaginary freq.
        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._max_iteration = max_iteration
        self._min_iteration = min_iteration
        self._traverse = traverse
        self._is_cell_relaxed = is_cell_relaxed

        self._stage = 0
        self._tasks = []

        self._energy = None
        self._cell = None
        self._phonon = None # Phonopy object
        self._phonon_fc3 = None # Phono3py object
        self._phonon_fc3_tasks = None

    def get_phonon(self):
        return self._phonon

    def get_phonon_fc3(self):
        for i, task in enumerate(self._phonon_fc3_tasks[1:]):
            forces_fc3.append(task.get_properties()['forces'][-1])
        disp_dataset = self._phonon_fc3.get_displacement_dataset()
        self._phonon_fc3.produce_fc3(forces_fc3)

        return self._phonon_fc3

    def get_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._phonon_fc3_tasks[0].get_cell()

    def get_energy(self):
        """Return energies at geometry optimization steps"""
        return self._energy

    def set_status(self):
        if self._stage == 0:
            task = self._tasks[0]
            if task.done():
                status = task.get_status()
                if status == "done":
                    self._status = "next"
                else:
                    self._status = status
        else:
            done = True
            terminate = False
            for i, task in enumerate(self._tasks):
                done &= task.done()
                terminate |= (task.get_status() == "terminate")

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
            self._phonon_fc3_tasks = [None]
            self._set_stage1()
        else:
            self._set_stage0()

    def end(self):
        pass

    def done(self):
        return (self._status == "terminate" or
                self._status == "done" or
                self._status == "max_iteration" or
                self._status == "next" or
                self._status == "imaginary_mode")

    def __next__(self):
        return self.next()

    def next(self):
        if self._stage == 0:
            if "next" in self._status:
                self._energy = self._tasks[0].get_energy()
                self._comment = "%s\\n%f" % (
                    self._tasks[0].get_space_group()['international'],
                    self._energy)
                self._set_stage1()
                return self._tasks
            elif "terminate" in self._status and self._traverse == "restart":
                self._traverse = False
                self._set_stage0()
                return self._tasks
            else:
                raise StopIteration
        elif self._stage == 1:
            if "next" in self._status:
                disp_dataset = self._phonon_fc3.get_displacement_dataset()
                for disp1, task in zip(disp_dataset['first_atoms'], self._tasks):
                    disp1['forces'] = task.get_properties()['forces'][-1]
                write_FORCE_SETS(disp_dataset)
                self._phonon.set_displacement_dataset(disp_dataset)
                self._phonon.produce_force_constants(
                    calculate_full_force_constants=False)
                if self._exist_imaginary_mode():
                    self._status = "imaginary_mode"
                    self._write_yaml()
                    self._tasks = []
                    raise StopIteration
                else:
                    self._set_stage2()
                    return self._tasks
            elif "terminate" in self._status and self._traverse == "restart":
                self._reset_stage1()
                return self._tasks
            else:
                raise StopIteration
        elif self._stage == 2:
            if "next" in self._status:
                self._status = "done"
                forces_fc3 = []
                for i, task in enumerate(self._phonon_fc3_tasks[1:]):
                    forces_fc3.append(task.get_properties()['forces'][-1])
                disp_dataset = self._phonon_fc3.get_displacement_dataset()
                write_FORCES_FC3(disp_dataset, forces_fc3)
                self._tasks = []
                raise StopIteration
            elif "terminate" in self._status and self._traverse == "restart":
                self._reset_stage2()
                return self._tasks
            else:
                raise StopIteration
        else: # stage2
            pass

    def _set_stage0(self):
        self._status = "equilibrium"
        task = self._get_equilibrium_task()
        self._phonon_fc3_tasks = [task]
        self._tasks = [task]

    def _set_stage1(self):
        self._set_phonon_fc3()
        if self._check_imaginary:
            self._stage = 1
            self._status = "fc2_displacements"
            disp_dataset = self._phonon_fc3.get_displacement_dataset()
            self._tasks = self._get_displacement_tasks(
                stop=len(disp_dataset['first_atoms']))
            self._phonon_fc3_tasks += self._tasks
        else:
            self._set_stage2()

    def _reset_stage1(self):
        self._traverse = False
        disp_terminated = []
        for i, task in enumerate(self._tasks):
            if task.get_status() == "terminate":
                disp_terminated.append(i)
        disp_dataset = self._phonon_fc3.get_displacement_dataset()
        tasks = self._get_displacement_tasks(
            stop=len(disp_dataset['first_atoms']))
        self._tasks = []
        for i in disp_terminated:
            self._tasks.append(tasks[i])
            self._phonon_fc3_tasks[i + 1] = tasks[i]
        self._status = "fc2_displacements"

    def _set_stage2(self):
        self._stage = 2
        self._status = "fc3_displacements"
        if self._check_imaginary:
            disp_dataset = self._phonon_fc3.get_displacement_dataset()
            start_index = len(disp_dataset['first_atoms'])
        else:
            start_index = 0
        self._tasks = self._get_displacement_tasks(start=start_index)
        self._phonon_fc3_tasks += self._tasks

    def _reset_stage2(self):
        self._traverse = False
        disp_terminated = []
        for i, task in enumerate(self._tasks):
            if task.get_status() == "terminate":
                disp_terminated.append(i)

        if self._check_imaginary:
            disp_dataset = self._phonon_fc3.get_displacement_dataset()
            start_index = len(disp_dataset['first_atoms'])
        else:
            start_index = 0
        tasks = self._get_displacement_tasks(start=start_index)
        self._tasks = []
        for i in disp_terminated:
            self._tasks.append(tasks[i])
            self._phonon_fc3_tasks[i + 1 + start_index] = tasks[i]
        self._status = "fc3_displacements"

    def _set_phonon_fc3(self):
        cell = self.get_cell()
        phonopy_cell = cell2atoms(cell)
        self._phonon = Phonopy(phonopy_cell,
                               self._supercell_matrix,
                               primitive_matrix=self._primitive_matrix,
                               dynamical_matrix_decimals=14,
                               force_constants_decimals=14)
        self._phonon_fc3 = Phono3py(phonopy_cell,
                                    self._supercell_matrix,
                                    primitive_matrix=self._primitive_matrix)
        self._phonon_fc3.generate_displacements(distance=self._distance,
                                                is_diagonal=self._is_diagonal)
        supercell = self._phonon_fc3.get_supercell()
        disp_dataset = self._phonon_fc3.get_displacement_dataset()
        self._phonon.set_displacement_dataset(disp_dataset)
        write_poscar(cell, "POSCAR-unitcell")
        write_disp_yaml(self._phonon.get_displacements(), supercell)
        write_disp_fc3_yaml(disp_dataset, supercell)

    def _exist_imaginary_mode(self):
        if self._primitive_matrix is None:
            pmat = np.eye(3)
        else:
            pmat = self._primitive_matrix
        exact_point_matrix = np.dot(np.linalg.inv(self._supercell_matrix),
                                    pmat).T
        max_integer = np.rint(np.amax(np.abs(np.linalg.inv(exact_point_matrix))))
        q_points = []
        for i in np.arange(-max_integer, max_integer + 1):
            for j in np.arange(-max_integer, max_integer + 1):
                for k in np.arange(-max_integer, max_integer + 1):
                    q = np.dot(exact_point_matrix, [i, j, k])
                    if (-1 < q).all() and (q < 1).all():
                        q_points.append(q)
        self._phonon.set_qpoints_phonon(q_points)
        frequencies = self._phonon.get_qpoints_phonon()[0]
        if (frequencies < self._cutoff_frequency).any():
            self._log = "Stop at phonon calculation due to imaginary modes"
            return True
        else:
            return False

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        if self._lattice_tolerance is not None:
            w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
        if self._stress_tolerance is not None:
            w.write("stress_tolerance: %f\n" % self._stress_tolerance)
            w.write("pressure_target: %f\n" % self._pressure_target)
        w.write("force_tolerance: %f\n" % self._force_tolerance)
        if self._max_increase is None:
            w.write("max_increase: unset\n")
        else:
            w.write("max_increase: %f\n" % self._max_increase)
        w.write("max_iteration: %d\n" % self._max_iteration)
        w.write("min_iteration: %d\n" % self._min_iteration)
        w.write("supercell_matrix:\n")
        for row in self._supercell_matrix:
            w.write("- [ %3d, %3d, %3d ]\n" % tuple(row))
        if self._primitive_matrix is not None:
            w.write("primitive_matrix:\n")
            for row in self._primitive_matrix:
                w.write("- [ %6.3f, %6.3f, %6.3f ]\n" % tuple(row))
        w.write("distance: %f\n" % self._distance)
        if self._phonon_fc3_tasks[0] is not None:
            w.write("iteration: %d\n" % self._phonon_fc3_tasks[0].get_stage())
            if self._energy:
                w.write("electric_total_energy: %20.10f\n" % self._energy)
        w.write("status: %s\n" % self._status)
        w.write("tasks:\n")
        for task in self._phonon_fc3_tasks:
            if task and task.get_status():
                w.write("- name:   %s\n" % task.get_name())
                w.write("  status: %s\n" % task.get_status())
        w.close()
