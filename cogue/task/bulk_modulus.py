from phonopy.qha import BulkModulus as PhonopyBulkModulus
from phonopy.units import EVAngstromToGPa

from cogue.crystal.cell import get_strained_cells
from cogue.task import TaskElement
from cogue.task.structure_optimization import StructureOptimizationYaml


class BulkModulusBase(TaskElement, StructureOptimizationYaml):
    """BulkModulus class

    Three stages:
    1. Structure optimization of input cell
    2. Total energiy calculations with +1 and -1 % volumes
    3. Calculate bulk modulus from stress of two cells
    or
    1. Structure optimization of input cell
    2. Total energy calculations at strains specified
    3. Calculate bulk modulus by fitting EOS

    """

    def __init__(
        self,
        directory=None,
        name=None,
        strains=None,
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
        self._task_type = "bulk_modulus"

        self._strains = strains

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
        self._bulk_modulus = None
        self._all_tasks = None

        self._eos = None  # Energy function of V, only for strains mode

    def get_bulk_modulus(self):
        return self._bulk_modulus

    def get_equation_of_state(self):
        return self._eos

    def get_cell(self):
        if self._is_cell_relaxed:
            return self._cell
        else:
            return self._all_tasks[0].get_cell()

    def get_all_tasks(self):
        return self._all_tasks

    def set_status(self):
        done = True
        terminate = False

        if self._stage == 0:
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

        self._write_yaml()

    def begin(self):
        if not self._job:
            print("set_job has to be executed.")
            raise RuntimeError

        if self._is_cell_relaxed:
            self._all_tasks = [None]
            self._prepare_next()
        else:
            self._set_stage0()

    def done(self):
        return (
            self._status == "done"
            or self._status == "terminate"
            or self._status == "max_iteration"
            or self._status == "next"
        )

    def __next__(self):
        return self.next()

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next()
                return self._tasks
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                self._set_stage0()
                return self._tasks
        else:
            if self._status == "next":
                self._status = "done"
                self._tasks = []
                if self._strains is None:
                    self._calculate_bulk_modulus_from_plus_minus()
                else:
                    self._calculate_bulk_modulus_by_fit_to_vinet()
            elif self._status == "terminate" and self._traverse == "restart":
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)
                cell = self.get_cell()
                if self._strains is None:
                    tasks = self._get_plus_minus_tasks(cell)
                else:
                    tasks = self._get_strained_cell_tasks(cell)
                self._tasks = []
                for i in terminated:
                    self._tasks.append(tasks[i])
                    self._all_tasks[i + 1] = tasks[i]
                self._status = "strains"
                return self._tasks

        self._write_yaml()
        raise StopIteration

    def _calculate_bulk_modulus_from_plus_minus(self):
        V = self.get_cell().get_volume()
        V_p = self._all_tasks[1].get_cell().get_volume()
        V_m = self._all_tasks[2].get_cell().get_volume()
        s_p = self._all_tasks[1].get_stress()
        s_m = self._all_tasks[2].get_stress()

        self._bulk_modulus = -(s_p - s_m).trace() / 3 * V / (V_p - V_m)

    def _calculate_bulk_modulus_by_fit_to_vinet(self):
        energies = [task.get_energy() for task in self._all_tasks[1:]]
        volumes = [task.get_cell().get_volume() for task in self._all_tasks[1:]]

        phonopy_bulk_modulus = PhonopyBulkModulus(volumes, energies)
        self._bulk_modulus = phonopy_bulk_modulus.get_bulk_modulus()
        self._bulk_modulus *= EVAngstromToGPa
        params = phonopy_bulk_modulus.get_parameters()

        def eos(v):
            return phonopy_bulk_modulus.get_eos()(v, *params)

        self._eos = eos

        with open("e-v.dat", "w") as w:
            w.write("#   cell volume        energy of cell " "other than phonon\n")
            for e, v in zip(energies, volumes):
                w.write("%20.13f %20.13f\n" % (v, e))

    def _set_stage0(self):
        self._stage = 0
        self._status = "equilibrium"
        task = self._get_equilibrium_task()
        self._all_tasks = [task]
        self._tasks = [task]

    def _prepare_next(self):
        cell = self.get_cell()

        self._stage = 1
        self._status = "strains"

        if self._strains is None:
            tasks = self._get_plus_minus_tasks(cell)
        else:
            tasks = self._get_strained_cell_tasks(cell)

        self._all_tasks += tasks
        self._tasks = tasks

    def _get_plus_minus_tasks(self, cell):
        cell_plus, cell_minus = get_strained_cells(cell, [0.01, -0.01])
        plus = self._get_equilibrium_task(
            index=1, cell=cell_plus, max_iteration=3, min_iteration=1, directory="plus"
        )
        minus = self._get_equilibrium_task(
            index=1,
            cell=cell_minus,
            max_iteration=3,
            min_iteration=1,
            directory="minus",
        )
        return plus, minus

    def _get_strained_cell_tasks(self, cell_orig):
        tasks = []
        for i, cell in enumerate(get_strained_cells(cell_orig, self._strains)):
            tasks.append(
                self._get_equilibrium_task(
                    index=1,
                    cell=cell,
                    max_iteration=3,
                    min_iteration=1,
                    directory="strain-%02d" % i,
                )
            )
        return tasks

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        lines += self._get_structopt_yaml_lines()
        cell = self.get_cell()
        if cell:
            lines += cell.get_yaml_lines()

        if self._bulk_modulus:
            lines.append("bulk_modulus: %f\n" % self._bulk_modulus)

        return lines
