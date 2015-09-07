from cogue.task import TaskElement

class BulkModulusBase(TaskElement):
    """BulkModulus class

    Three stages:
    1. structure optimization of input cell
    2. create cells with +1% and -1% volume and optimize them
    3. calculate bulk modulus from stress of two cells
    
    """
    
    def __init__(self,
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
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "bulk_modulus"
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
        self._bm_tasks = None

    def get_bulk_modulus(self):
        return self._bulk_modulus

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
        
    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            self._bm_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._bm_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._bm_tasks[0]]

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "max_iteration" or
                self._status == "next")

    def next(self):    
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._bm_tasks[0].get_cell())
                return self._tasks
        else:
            if self._status == "next":
                stress_p = self._bm_tasks[1].get_stress()
                stress_m = self._bm_tasks[2].get_stress()

                if (stress_p is None or stress_m is None):
                    self._status = "terminate"
                else:
                    self._calculate_bulk_modulus()
                    self._status = "done"

        self._write_yaml()
        raise StopIteration

    def _calculate_bulk_modulus(self):
        if self._is_cell_relaxed:
            V = self._cell.get_volume()
        else:
            V = self._bm_tasks[0].get_cell().get_volume()
        V_p = self._bm_tasks[1].get_cell().get_volume()
        V_m = self._bm_tasks[2].get_cell().get_volume()
        s_p = self._bm_tasks[1].get_stress()
        s_m = self._bm_tasks[2].get_stress()

        self._bulk_modulus = - (s_p - s_m).trace() / 3 * V / (V_p - V_m)
        
    def _prepare_next(self, cell):
        self._stage = 1
        self._status = "plus minus"
        plus, minus = self._get_plus_minus_tasks(cell)
        self._bm_tasks.append(plus)
        self._bm_tasks.append(minus)
        self._tasks = self._bm_tasks[1:]

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
        if self._bm_tasks[0] is not None:
            w.write("iteration: %d\n" % self._bm_tasks[0].get_stage())
        w.write("status: %s\n" % self._status)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._bm_tasks[0].get_cell()

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

        if self._bulk_modulus:
            w.write("bulk_modulus: %f\n" % self._bulk_modulus)
        w.close()
