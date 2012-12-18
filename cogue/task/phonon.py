from cogue.task import TaskElement
from cogue.crystal.converter import atoms2cell, write_v_sim
from phonopy import Phonopy
from phonopy.structure.atoms import Atoms
from phonopy.hphonopy.file_IO import write_disp_yaml
from phonopy.hphonopy.file_IO import write_FORCE_SETS
from cogue.crystal.vasp_io import write_poscar

class PhononBase(TaskElement):
    """PhononBase class

    This is an interface to phonopy.

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
        self._task_type = "phonon"
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

    def get_phonon(self):
        return self._phonon

    def get_equilibrium_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._phonon_tasks[0].get_cell()

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

        self._overwrite_settings()

        if self._is_cell_relaxed:
            self._phonon_tasks = [None]
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
        else: # task 1..n: displaced supercells
            if "next" in self._status:
                self._status = "done"
                forces = []
                for task in self._phonon_tasks[1:]:
                    forces.append(task.get_properties()['forces'][-1])
                self._write_FORCE_SETS(forces)
                self._phonon.set_post_process(self._primitive_matrix,
                                              forces,
                                              force_constants_decimals=14)
                self._tasks = []
                raise StopIteration
            elif "terminate" in self._status and self._traverse == "restart":
                self._traverse = False
                disp_terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        disp_terminated.append(i)
                tasks = self._get_displacement_tasks()[1:]
                self._tasks = []
                for i in disp_terminated:
                    self._tasks.append(tasks[i])
                    self._phonon_tasks[i + 1] = tasks[i]
                self._status = "displacements"
                return self._tasks
            else:
                raise StopIteration

    def _set_stage0(self):
        self._status = "equilibrium"
        task = self._get_equilibrium_task()
        self._phonon_tasks = [task]
        self._tasks = [task]
        
    def _set_stage1(self):
        self._stage = 1
        self._status = "displacements"
        self._set_phonon()
        self._tasks = self._get_displacement_tasks()[1:]
        self._phonon_tasks += self._tasks

    def _set_phonon(self):
        cell = self.get_equilibrium_cell()
        phonopy_cell = Atoms(
            cell=cell.get_lattice().T,
            scaled_positions=cell.get_points().T,
            symbols=cell.get_symbols())
        
        self._phonon = Phonopy(phonopy_cell,
                               self._supercell_matrix,
                               is_auto_displacements=False)
        self._phonon.generate_displacements(distance=self._distance,
                                            is_diagonal=False)

        supercell = self._phonon.get_supercell()
        displacements = self._phonon.get_displacements()

        write_poscar(cell, "POSCAR-unitcell")
        write_disp_yaml(displacements, supercell)

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
            w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
            w.write("pressure_target: %f\n" % self._pressure_target)
            w.write("stress_tolerance: %f\n" % self._stress_tolerance)
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


