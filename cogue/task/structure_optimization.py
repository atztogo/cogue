import os
from cogue.task import TaskElement
from cogue.crystal.symmetry import \
    get_crystallographic_cell, get_symmetry_dataset
from cogue.crystal.converter import get_primitive

class StructureOptimizationBase(TaskElement):
    def __init__(self,
                 directory="structure_optimization",
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
        self._task_type = "structure_optimization"

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
        self._so_tasks = None
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

        self._so_tasks = [task]
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
            if not stress == None:
                self._stress = stress
            if not forces == None:
                self._forces = forces
            if not energy == None:
                self._energy = energy

        if self._status == "terminate" and self._traverse == "restart":
            self._traverse = False
            if self._stage > 2:
                self._stage -= 2
                task = self._so_tasks.pop()
                task = self._so_tasks.pop()
                self._next_cell = task.get_cell()
            else:
                self._so_tasks = []
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

    def _set_next_task(self):
        self._stage += 1
        self._status = "stage %d" % self._stage
        task = self._get_next_task(self._next_cell)
        self._so_tasks.append(task)
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
        w.write("iteration: %d\n" % self._stage)
        w.write("status: %s\n" % self._status)

        cell = self._so_tasks[-1].get_current_cell()
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

            if self._space_group:
                w.write("symmetry_tolerance: %s\n" % self._symmetry_tolerance)
                w.write("space_group_type: %s\n" %
                        self._space_group['international_standard'])
                w.write("space_group_number: %d\n" %
                        self._space_group['number'])

        if not self._energy == None:
            w.write("energy: %20.10f\n" % self._energy)

        if not self._forces == None:
            w.write("forces:\n")
            for i, v in enumerate(self._forces):
                w.write("- [ %15.10f, %15.10f, %15.10f ] # %d\n" %
                        (v[0], v[1], v[2], i + 1))

        if not self._stress == None:
            w.write("stress:\n")
            for x, v in zip(('x', 'y', 'z'), self._stress):
                w.write("- [ %15.10f, %15.10f, %15.10f ] # %sx %sy %sz\n" % (v[0], v[1], v[2], x, x, x))

        w.close()
        
