from cogue.task import TaskElement
from cogue.task.phonon import PhononYaml
import numpy as np

class ModeGruneisenBase(TaskElement, PhononYaml):
    """ModeGruneisen class

    Three stages:
    1. structure optimization of input cell
    2. create cells with, e.g. +1% and -1% volumes or strains, and optimize them
    3. calculate phonons for three cells
    4. calculate mode-Gruneisen parameters (This is not yet implemented.)
    
    """
    
    def __init__(self,
                 directory=None,
                 name=None,
                 delta_strain=None,
                 strain=None,
                 bias=None,
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
                 is_cell_relaxed=False,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "mode_gruneisen"

        self._bias = bias
        if delta_strain is None:
            self._delta_strain = 0.001
        else:
            self._delta_strain = delta_strain
        (self._delta_strain_minus,
         self._delta_strain_orig,
         self._delta_strain_plus) = self._get_delta_strains()
        if strain is None:
            self._strain = np.eye(3)
        else:
            self._strain = self._get_strain(strain)
        
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
        self._tasks = None

        self._cell = None
        self._mode_gruneisen = None
        self._all_tasks = None

    def get_mode_gruneisen(self):
        return self._mode_gruneisen

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
            self._all_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._all_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._all_tasks[0]]

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):    
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._all_tasks[0].get_cell())
                return self._tasks
        else:
            if self._status == "next":
                self._calculate_mode_gruneisen()
                self._status = "done"
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)
                tasks = self._get_phonon_tasks(cell)
                self._tasks = []
                for i in terminated:
                    self._tasks.append(tasks[i])
                    self._all_tasks[i + 1] = tasks[i]
                self._status = "phonons"
                return self._tasks

        self._write_yaml()
        raise StopIteration

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._all_tasks[0].get_cell()
        lines += self._get_phonon_yaml_lines(cell)

        return lines

    def _calculate_mode_gruneisen(self):
        self._mode_gruneisen = None
        
    def _prepare_next(self, cell):
        self._stage = 1
        self._status = "phonons"
        tasks = self._get_phonon_tasks(cell)
        self._tasks = tasks
        self._all_tasks += tasks

    def _get_delta_strains(self):
        if self._bias is "plus":
            return (np.eye(3),
                    self._get_strain(self._delta_strain, factor=1),
                    self._get_strain(self._delta_strain, factor=2))
        elif self._bias is "minus":
            return (self._get_strain(self._delta_strain, factor=-2),
                    self._get_strain(self._delta_strain, factor=-1),       
                    np.eye(3))
        else:
            return (self._get_strain(self._delta_strain, factor=-1),
                    np.eye(3),
                    self._get_strain(self._delta_strain, factor=1))
        
    def _get_strain(self, strain, factor=1):
        if isinstance(strain, int) or isinstance(strain, float):
            return (1 + factor * strain) ** (1.0 / 3) * np.eye(3)
        else:
            return np.eye(3) + factor * np.array(strain)

    def _get_phonon_tasks(self, cell):
        lattice = np.dot(self._strain, cell.get_lattice())
        cell_orig = cell.copy()
        cell_orig.set_lattice(np.dot(self._delta_strain_orig, lattice))
        cell_minus = cell.copy()
        cell_minus.set_lattice(np.dot(self._delta_strain_minus, lattice))
        cell_plus = cell.copy()
        cell_plus.set_lattice(np.dot(self._delta_strain_plus, lattice))

        minus = self._get_phonon_task(cell_minus,
                                      "minus",
                                      is_cell_relaxed=(
                                          self._is_cell_relaxed and
                                          self._bias == "plus"))

        orig = self._get_phonon_task(cell_orig,
                                     "orig",
                                     is_cell_relaxed=(
                                         self._is_cell_relaxed and
                                         self._bias != "plus" and
                                         self._bias != "minus"))

        plus = self._get_phonon_task(cell_plus,
                                     "plus",
                                     is_cell_relaxed=(
                                         self._is_cell_relaxed and
                                         self._bias == "minus"))

        return orig, plus, minus
