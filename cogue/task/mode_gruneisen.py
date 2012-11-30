from cogue.task import TaskElement

class ModeGruneisenBase(TaskElement):
    """ModeGruneisen class

    Three stages:
    1. structure optimization of input cell
    2. create cells with +1% and -1% volume and optimize them
    3. calculate phonons for three cells
    4. calculate mode-Gruneisen parameters
    
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
        self._task_type = "mode_gruneisen"

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
        self._mg_tasks = None

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

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            self._mg_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._mg_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._mg_tasks[0]]

    def end(self):
        self._write_yaml()

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):    
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._mg_tasks[0].get_cell())
            else:
                raise StopIteration
        else:
            if self._status == "next":
                self._calculate_mode_gruneisen()
                self._status = "done"
            raise StopIteration

        return self._tasks

    def _calculate_mode_gruneisen(self):
        self._mode_gruneisen = None
        
    def _prepare_next(self, cell):
        self._stage = 1
        self._status = "phonons"
        self._mg_tasks += self._get_phonon_tasks(cell)
        self._tasks = self._mg_tasks[1:]

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        if self._mg_tasks[0]:
            w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
            w.write("pressure_target: %f\n" % self._pressure_target)
            w.write("stress_tolerance: %f\n" % self._stress_tolerance)
            w.write("force_tolerance: %f\n" % self._force_tolerance)
            w.write("max_increase: %f\n" % self._max_increase)
            w.write("max_iteration: %d\n" % self._max_iteration)
            w.write("min_iteration: %d\n" % self._min_iteration)
            w.write("iteration: %d\n" % self._mg_tasks[0].get_stage())
        w.write("status: %s\n" % self._status)
        w.write("supercell_matrix:\n")
        for row in self._supercell_matrix:
            w.write("- [ %3d, %3d, %3d ]\n" % tuple(row))
        w.write("primitive_matrix:\n")
        for row in self._primitive_matrix:
            w.write("- [ %6.3f, %6.3f, %6.3f ]\n" % tuple(row))
        w.write("distance: %f\n" % self._distance)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._mg_tasks[0].get_cell()

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
