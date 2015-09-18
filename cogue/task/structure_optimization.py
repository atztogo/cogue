import os
from cogue.task import TaskElement
from cogue.task.oneshot_calculation import OneShotCalculationYaml
from cogue.crystal.symmetry import get_crystallographic_cell, get_symmetry_dataset
from cogue.crystal.converter import get_primitive

class StructureOptimizationYaml(OneShotCalculationYaml):
    def _get_structopt_yaml_lines(self):
        lines = []
        if self._all_tasks:
            lines.append("tasks:")
        for task in self._all_tasks:
            if not task:
                continue
            for i, line in enumerate(TaskElement.get_yaml_lines(task)):
                if i == 0:
                    lines.append("- " + line)
                else:
                    lines.append("  " + line)
        if self._lattice_tolerance is not None:
            lines.append("lattice_tolerance: %f" % self._lattice_tolerance)
        if self._stress_tolerance is not None:
            lines.append("stress_tolerance: %f" % self._stress_tolerance)
            lines.append("pressure_target: %f" % self._pressure_target)
        lines.append("force_tolerance: %f" % self._force_tolerance)
        if self._max_increase is None:
            lines.append("max_increase: unset")
        else:
            lines.append("max_increase: %f" % self._max_increase)
        lines.append("max_iteration: %d" % self._max_iteration)
        lines.append("min_iteration: %d" % self._min_iteration)

        return lines

class StructureOptimizationBase(TaskElement, StructureOptimizationYaml):
    def __init__(self,
                 directory="structopt",
                 name=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 max_iteration=None,
                 min_iteration=None,
                 find_symmetry=True,
                 impose_symmetry=False,
                 symmetry_tolerance=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "structopt"

        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._max_iteration = max_iteration
        self._min_iteration = min_iteration
        self._find_symmetry = find_symmetry
        self._impose_symmetry = impose_symmetry
        self._symmetry_tolerance = symmetry_tolerance
        self._traverse = traverse

        self._stage = 1
        self._tasks = None

        self._cell = None
        self._next_cell = None
        self._all_tasks = None
        self._stress = None
        self._forces = None
        self._energy = None
        self._space_group = None

    def get_symmetry_tolerance(self):
        return self._symmetry_tolerance

    def get_cell(self):
        return self._next_cell

    def get_initial_cell(self):
        return self._cell

    def get_stage(self):
        return self._stage

    def get_stress(self):
        return self._stress

    def get_forces(self):
        return self._forces

    def get_energy(self):
        return self._energy

    def get_space_group(self):
        return self._space_group

    def set_status(self):
        task = self._tasks[0]
        if task.done():
            self._status = task.get_status()

        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        self._overwrite_settings()

        self._status = "stage 1"
        if self._impose_symmetry:
            prim_cell = get_primitive(self._cell,
                                      tolerance=self._symmetry_tolerance)
            self._space_group = get_symmetry_dataset(prim_cell)
            task = self._get_next_task(prim_cell)
        else:
            if self._find_symmetry:
                self._space_group = get_symmetry_dataset(
                    self._cell,
                    self._symmetry_tolerance)
            task = self._get_next_task(self._cell)
        if self._space_group:
            self._comment = self._space_group['international_standard']

        self._all_tasks = [task]
        self._tasks = [task]
        self._write_yaml()

    def done(self):
        return (self._status == "next" or
                self._status == "done" or
                self._status == "terminate" or
                self._status == "max_iteration")

    def next(self):
        if self._status == "terminate":
            self._stress = None
            self._forces = None
            self._energy = None
            self._next_cell = None
        else:
            task = self._tasks[0]
            self._next_cell = task.get_current_cell()
            stress = task.get_stress()
            forces = task.get_forces()
            energy = task.get_energy()
            if stress is not None:
                self._stress = stress
            if forces is not None:
                self._forces = forces
            if energy is not None:
                self._energy = energy

        if self._status == "terminate" and self._traverse == "restart":
            self._traverse = False
            if self._stage > 2:
                self._stage -= 2
                task = self._all_tasks.pop()
                task = self._all_tasks.pop()
                self._next_cell = task.get_cell()
            else:
                self._all_tasks = []
                self._stage = 0
                self._next_cell = self._cell
            self._status = "next"
                
        if self._next_cell:
            if self._impose_symmetry:
                self._next_cell = get_primitive(
                    self._next_cell,
                    tolerance=self._symmetry_tolerance)
                self._space_group = get_symmetry_dataset(self._next_cell)
            else:
                if self._find_symmetry:
                    self._space_group = get_symmetry_dataset(
                        self._next_cell,
                        tolerance=self._symmetry_tolerance)
            if self._space_group:
                self._comment = self._space_group['international_standard']
            
        if self._status == "done":
            if self._stage < self._min_iteration:
                self._status = "next"

        if self._status == "next":
            if self._stage == self._max_iteration:
                self._status = "max_iteration"
            else:
                self._set_next_task()

        self._write_yaml()
        if "stage" in self._status:
            return self._tasks
        else:
            self._tasks = []
            raise StopIteration

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        lines.append("iteration: %d" % self._stage)
        lines += self._get_structopt_yaml_lines()
        cell = self._all_tasks[-1].get_current_cell()
        lines += self._get_oneshot_yaml_lines(cell)
        if self._space_group:
            lines.append("symmetry_tolerance: %s" %
                         self._symmetry_tolerance)
            lines.append("space_group_type: %s" %
                         self._space_group['international_standard'])
            lines.append("space_group_number: %d" %
                         self._space_group['number'])

        return lines

    def _set_next_task(self):
        self._stage += 1
        self._status = "stage %d" % self._stage
        task = self._get_next_task(self._next_cell)
        self._all_tasks.append(task)
        self._tasks = [task]
