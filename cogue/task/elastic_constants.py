from cogue.task import TaskElement
from cogue.task.structure_optimization import StructureOptimizationYaml


class ElasticConstantsYaml(StructureOptimizationYaml):
    def _get_ec_yaml_lines(self, cell):
        lines = self._get_structopt_yaml_lines()
        if self._elastic_constants is not None:
            lines.append("elastic_constants:")
            for v in self._elastic_constants:
                lines.append(
                    "- [ %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f ]" % tuple(v)
                )
        if cell:
            lines += cell.get_yaml_lines()

        return lines


class ElasticConstantsBase(TaskElement, ElasticConstantsYaml):
    """ElasticConstantsBase class

    Calculate elastic constants
    1. Structure optimization to obtain equilibrium structure
    2. Calculate elastic constants by a task that returns elastic constants.

    """

    def __init__(
        self,
        directory=None,
        name=None,
        lattice_tolerance=None,
        force_tolerance=None,
        pressure_target=None,
        stress_tolerance=None,
        max_increase=None,
        max_iteration=None,
        min_iteration=None,
        is_cell_relaxed=False,
        traverse=False,
    ):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "elastic_constants"
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
        self._all_tasks = None
        self._elastic_constants = None
        self._energy = None

    def get_elastic_constants(self):
        return self._elastic_constants

    def get_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._all_tasks[0].get_cell()

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
        return (
            self._status == "done"
            or self._status == "terminate"
            or self._status == "next"
        )

    def __next__(self):
        return self.next()

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._energy = self._tasks[0].get_energy()
                self._set_stage1()
                return self._tasks
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                self._set_stage0()
                return self._tasks
        else:
            if self._status == "next":
                self._elastic_constants = self._all_tasks[1].get_elastic_constants()
                self._status = "done"
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                self._set_stage1()
                return self._tasks

        self._write_yaml()
        raise StopIteration

    def _set_stage0(self):
        self._status = "equilibrium"
        task = self._get_equilibrium_task()
        self._all_tasks = [task]
        self._tasks = [task]

    def _set_stage1(self):
        self._stage = 1
        self._status = "elastic constants"
        ec_task = self._get_ec_task(self.get_cell())
        if len(self._all_tasks) == 2:
            self._all_tasks[1] = ec_task
        else:
            self._all_tasks.append(ec_task)
        self._tasks = [self._all_tasks[1]]

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        cell = self.get_cell()
        lines += self._get_ec_yaml_lines(cell)
        if self._all_tasks[0] is not None:
            if self._energy:
                lines.append("electric_total_energy: %20.10f" % self._energy)

        return lines
