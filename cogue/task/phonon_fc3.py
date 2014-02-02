from cogue.task import TaskElement
from anharmonic.phonon3 import Phono3py
from phonopy.structure.atoms import Atoms
from anharmonic.file_IO import write_disp_fc3_yaml
from anharmonic.file_IO import write_FORCES_FC3
from cogue.interface.vasp_io import write_poscar
from cogue.task.phonon import cell2atoms

class PhononFC3Base(TaskElement):
    """PhononFC3Base class

    This is an interface to anharmonic phonopy.

    """

    def __init__(self,
                 directory=None,
                 name=None,
                 supercell_matrix=None,
                 primitive_matrix=None,
                 distance=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 traverse=False,
                 is_cell_relaxed=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "anharmonic_phonon"
        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._distance = distance
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
        done = True
        terminate = False
        for task in self._tasks:
            done &= task.done()
            if task.get_status() == "terminate":
                terminate = True
        if done:
            if terminate:
                self._status = "terminate"
            else:
                self._status = "next"

        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            self._phonon_fc3_tasks = [None]
            self._set_stage1()
        else:
            self._set_stage0()

    def end(self):
        pass

    def done(self):
        return ("terminate" in self._status or 
                "done" in self._status or
                "next" in self._status)

    def next(self):
        if self._stage == 0:
            if "next" in self._status:
                self._energy = self._tasks[0].get_energy()
                self._comment = "%s\\n%f" % (
                    self._tasks[0].get_space_group()['international_standard'],
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
                self._status = "done"
                disp_dataset = self._phonon_fc3.get_displacements()
                for disp1, task in zip(disp_dataset['first_atom'], self._tasks):
                    disp1['forces'] = task.get_properties()['forces'][-1]
                self._phonon.set_displacement_dataset(disp_dataset)
                write_FORCE_SETS_from_dataset(disp_dataset)
                self._phonon.produce_force_constants()
                self._set_stage2()
                return self._tasks
            elif "terminate" in self._status and self._traverse == "restart":
                self._traverse = False
                disp_terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        disp_terminated.append(i)
                disp_dataset = self._phonon_fc3.get_displacements()
                num_fc2_displacements = len(disp_dataset['first_atoms'])
                tasks = self._get_displacement_tasks(stop=num_fc2_displacements)
                self._tasks = []
                for i in disp_terminated:
                    self._tasks.append(tasks[i])
                    self._phonon_fc3_tasks[i + 1] = tasks[i]
                self._status = "fc2_displacements"
                return self._tasks
            else:
                raise StopIteration
        elif self._stage == 2:
            if "next" in self._status:
                self._status = "done"
                forces = []
                for task in self._tasks:
                    forces.append(task.get_properties()['forces'][-1])
                self._tasks = []
                raise StopIteration
            elif "terminate" in self._status and self._traverse == "restart":
                self._traverse = False
                disp_terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        disp_terminated.append(i)
                disp_dataset = self._phonon_fc3.get_displacements()
                num_fc2_displacements = len(disp_dataset['first_atoms'])
                tasks = self._get_displacement_tasks(stop=num_fc2_displacements)
                self._tasks = []
                for i in disp_terminated:
                    self._tasks.append(tasks[i])
                    self._phonon_fc3_tasks[i + 1] = tasks[i]
                self._status = "fc3_displacements"
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
        self._stage = 1
        self._status = "fc2_displacements"
        self._set_phonon_fc3()
        disp_dataset = self._phonon_fc3.get_displacements()
        num_fc2_displacements = len(disp_dataset['first_atoms'])
        self._tasks = self._get_displacement_tasks(stop=num_fc2_displacements)
        self._phonon_fc3_tasks += self._tasks

    def _set_stage2(self):
        self._stage = 2
        self._status = "fc3_displacements"
        disp_dataset = self._phonon_fc3.get_displacements()
        num_fc2_displacements = len(disp_dataset['first_atoms'])
        self._tasks = self._get_displacement_tasks(start=num_fc2_displacements)
        self._phonon_fc3_tasks += self._tasks
        
    def _set_phonon_fc3(self):
        cell = self.get_cell()
        phonopy_cell = cell2atoms(cell)
        self._phonon = Phonopy(phonopy_cell,
                               self._supercell_matrix,
                               primitive_matrix=self._primitive_matrix,
                               is_auto_displacements=False,
                               dynamical_matrix_decimals=14,
                               force_constants_decimals=14)
        self._phonon_fc3 = Phono3py(phonopy_cell,
                                    self._supercell_matrix,
                                    primitive_matrix=self._primitive_matrix)
        self._phonon_fc3.generate_displacements(distance=self._distance)
        supercell = self._phonon_fc3.get_supercell()
        disp_dataset = self._phonon_fc3.get_displacements()
        write_poscar(cell, "POSCAR-unitcell")
        write_disp_fc3_yaml(disp_dataset, supercell)

    def _write_FORCE_SETS(self, forces):
        displacements = [[x[0], x[1:4]]
                         for x in self._phonon.get_displacements()]
        natom = self._phonon.get_supercell().get_number_of_atoms()
        write_FORCE_SETS("FORCE_SETS",
                         natom,
                         displacements,
                         forces,
                         verbose=False)

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        if self._phonon_tasks[0]:
            if self._lattice_tolerance is not None:
                w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
            if self._stress_tolerance is not None:
                w.write("stress_tolerance: %f\n" % self._stress_tolerance)
                w.write("pressure_target: %f\n" % self._pressure_target)
            w.write("force_tolerance: %f\n" % self._force_tolerance)
            w.write("max_increase: %f\n" % self._max_increase)
            w.write("max_iteration: %d\n" % self._max_iteration)
            w.write("min_iteration: %d\n" % self._min_iteration)
            w.write("supercell_matrix:\n")
            for row in self._supercell_matrix:
                w.write("- [ %3d, %3d, %3d ]\n" % tuple(row))
            w.write("primitive_matrix:\n")
            for row in self._primitive_matrix:
                w.write("- [ %6.3f, %6.3f, %6.3f ]\n" % tuple(row))
            w.write("distance: %f\n" % self._distance)
            w.write("iteration: %d\n" % self._phonon_tasks[0].get_stage())
            if self._energy:
                w.write("electric_total_energy: %20.10f\n" % self._energy)
        w.write("status: %s\n" % self._status)
        w.write("tasks:\n")
        for task in self._phonon_tasks:
            if task and task.get_status():
                w.write("- name:   %s\n" % task.get_name())
                w.write("  status: %s\n" % task.get_status())
        w.close()


