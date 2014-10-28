from cogue.task import TaskElement
from phonopy import Phonopy
from phonopy.structure.atoms import Atoms
from phonopy.file_IO import write_disp_yaml
from phonopy.file_IO import write_FORCE_SETS
from cogue.interface.vasp_io import write_poscar
from cogue.crystal.cell import Cell

#########################
# Cell to Phonopy Atoms #
#########################
def cell2atoms(cell):
    return Atoms(cell=cell.get_lattice().T,
                 scaled_positions=cell.get_points().T,
                 masses=cell.get_masses(),
                 magmoms=cell.get_magnetic_moments(),
                 symbols=cell.get_symbols())

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
                 stop_condition=None,
                 traverse=False):

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
        self._stop_condition = stop_condition
        self._traverse = traverse

        self._stage = 0
        self._tasks = []

        self._energy = None
        self._space_group = None
        self._cell = None
        self._phonon = None # Phonopy object
        self._phonon_tasks = None # Phonopy object

    def get_phonon(self):
        return self._phonon

    def get_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._phonon_tasks[0].get_cell()

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
            terminate = True
            for task in self._tasks:
                done &= task.done()
                terminate &= (task.get_status() == "terminate")
                
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

    def done(self):
        return (self._status == "terminate" or
                self._status == "done" or
                self._status == "max_iteration" or
                self._status == "low_symmetry" or
                self._status == "next")

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._energy = self._tasks[0].get_energy()
                num_atom = len(self._tasks[0].get_cell().get_symbols())
                self._comment = self._space_group['international_standard']
                self._comment += "\\n%f/%d" % (self._energy, num_atom)
                self._set_stage1()
                return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._set_stage0()
                return self._tasks
                
        else: # task 1..n: displaced supercells
            if self._status == "next":
                self._status = "done"
                forces = []
                for task in self._tasks:
                    forces.append(task.get_properties()['forces'][-1])
                self._phonon.produce_force_constants(forces)
                write_FORCE_SETS(self._phonon.get_displacement_dataset())
                self._tasks = []
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
                    self._phonon_tasks[i + 1] = tasks[i]
                self._status = "displacements"
                return self._tasks
                
        self._write_yaml()
            
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
        self._tasks = self._get_displacement_tasks()
        self._phonon_tasks += self._tasks

    def _evaluate_stop_condition(self):
        if self._stop_condition:
            if "symmetry_operations" in self._stop_condition:
                num_ops = len(self._space_group['rotations'])
                if (num_ops <
                    self._stop_condition['symmetry_operations']):
                    self._status = "low_symmetry"
                    return True
                    
        return False
        
    def _set_phonon(self):
        cell = self.get_cell()
        phonopy_cell = cell2atoms(cell)
        self._phonon = Phonopy(phonopy_cell,
                               self._supercell_matrix,
                               primitive_matrix=self._primitive_matrix,
                               is_auto_displacements=False,
                               dynamical_matrix_decimals=14,
                               force_constants_decimals=14)
        self._phonon.generate_displacements(
            distance=self._distance,
            is_plusminus=self._displace_plusminus,
            is_diagonal=self._displace_diagonal)
        supercell = self._phonon.get_supercell()
        displacements = self._phonon.get_displacements()

        write_poscar(cell, "POSCAR-unitcell")
        write_disp_yaml(displacements, supercell)

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        w.write("supercell_matrix:\n")
        for row in self._supercell_matrix:
            w.write("- [ %3d, %3d, %3d ]\n" % tuple(row))
        w.write("primitive_matrix:\n")
        for row in self._primitive_matrix:
            w.write("- [ %6.3f, %6.3f, %6.3f ]\n" % tuple(row))
        w.write("distance: %f\n" % self._distance)
        
        if self._phonon_tasks[0]:
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


