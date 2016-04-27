import numpy as np
from cogue.task import TaskElement

class DensityOfStatesBase(TaskElement):
    """DensityOfStates class

    Three stages:
    1. structure optimization of input cell
    2. calculate charge density
    3. calculate density of states

    """

    def __init__(self,
                 directory=None,
                 name=None,
                 is_partial_dos=False,
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
        self._task_type = "density_of_states"
        self._is_partial_dos = is_partial_dos
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
        self._density_of_states = None
        self._dos_tasks = None

    def get_density_of_states(self):
        return self._density_of_states

    def set_status(self):
        done = True
        terminate = False

        if self._stage == 0 or self._stage == 1:
            task = self._tasks[0]
            if task.done():
                status = task.get_status()
                if status == "terminate" or status == "max_iteration":
                    self._status = status
                else:
                    self._status = "next"
        else:
            for task in self._tasks:
                done &= task.done()
                if task.get_status() == "terminate":
                    terminate = True
            if done:
                if terminate:
                    self._status = "terminate"
                else:
                    self._status = "next"

    def begin(self):
        if not self._job:
            print("set_job has to be executed.")
            raise

        if self._is_cell_relaxed:
            self._dos_tasks = [None]
            self._set_stage1()
        else:
            self._status = "equilibrium"
            self._dos_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._dos_tasks[0]]

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "max_iteration" or
                self._status == "next")

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._set_stage1()
                return self._tasks
        elif self._stage == 1:
            if self._status == "next":
                self._set_stage2()
                return self._tasks
        else:
            if self._status == "next":
                self._create_density_of_states()
                self._status = "done"
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                if len(self._dos_tasks) > 1:
                    self._dos_tasks = self._dos_tasks[:2]
                self._set_stage2()
                return self._tasks

        self._write_yaml()
        raise StopIteration

    def _create_density_of_states(self):
        pass

    def _set_stage1(self):
        if self._dos_tasks[0] is None:
            cell = self._cell
        else:
            cell = self._dos_tasks[0].get_cell()
        self._stage = 1
        self._status = "charge density"
        task = self._get_charge_density_task(cell)
        self._dos_tasks.append(task)
        self._tasks = [task]

    def _set_stage2(self):
        cell = self._dos_tasks[1].get_cell()
        properties = self._dos_tasks[1].get_properties()
        self._stage = 2
        self._status = "dos"
        task = self._get_dos_task(cell, properties=properties)
        self._dos_tasks.append(task)
        self._tasks = [task]

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
        if self._dos_tasks[0] is not None:
            w.write("iteration: %d\n" % self._dos_tasks[0].get_stage())
        w.write("status: %s\n" % self._status)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._dos_tasks[0].get_cell()

        if cell:
            lattice = cell.get_lattice().T
            points = cell.get_points().T
            symbols = cell.get_symbols()

            w.write("lattice:\n")
            for v, a in zip(lattice, ('a', 'b', 'c')) :
                w.write("- [ %22.16f, %22.16f, %22.16f ] # %s\n" %
                        (v[0], v[1], v[2], a))

            w.write("points:\n")
            for i, v in enumerate(points):
                w.write("- [ %20.16f, %20.16f, %20.16f ] # %d\n" %
                        (v[0], v[1], v[2], i + 1))

            w.write("symbols:\n")
            for i, v in enumerate(symbols):
                w.write("- %2s # %d\n" % (v, i + 1))
        w.close()
