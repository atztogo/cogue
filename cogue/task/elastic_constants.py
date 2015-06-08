from cogue.task import TaskElement

class ElasticConstantsBase(TaskElement):
    """ElasticConstantsBase class

    Calculate elastic constants
    1. Structure optimization to obtain equilibrium structure
    2. Calculate elastic constants by
    a task that returns elastic constants.

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
        self._ec_tasks = None
        self._elastic_constants = None

    def get_elastic_constants(self):
        return self._elastic_constants

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
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            self._ec_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._ec_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._ec_tasks[0]]

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._ec_tasks[0].get_cell())
                return self._tasks
        else:
            if self._status == "next":
                self._elastic_constants = \
                    self._ec_tasks[1].get_elastic_constants()
                self._status = "done"

        self._write_yaml()
        raise StopIteration

    def _prepare_next(self, cell):
        self._stage = 1
        self._status = "elastic constants"
        self._ec_tasks.append(self._get_ec_task(cell ))
        self._tasks = [self._ec_tasks[1]]

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        if self._ec_tasks[0]:
            if self._lattice_tolerance is not None:
                w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
            if self._stress_tolerance is not None:
                w.write("pressure_target: %f\n" % self._pressure_target)
                w.write("stress_tolerance: %f\n" % self._stress_tolerance)
            w.write("force_tolerance: %f\n" % self._force_tolerance)
            if self._max_increase is None:
                w.write("max_increase: unset\n")
            else:
                w.write("max_increase: %f\n" % self._max_increase)
            w.write("max_iteration: %d\n" % self._max_iteration)
            w.write("iteration: %d\n" % self._ec_tasks[0].get_stage())
        w.write("status: %s\n" % self._status)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._ec_tasks[0].get_cell()

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

        if self._elastic_constants is not None:
            w.write("elastic_constants:\n")
            for v in self._elastic_constants:
                w.write("- [ %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f ]\n" % tuple(v))
        w.close()

