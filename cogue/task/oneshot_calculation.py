import os
from cogue.task import TaskElement

class OneShotCalculationYaml:
    def _get_oneshot_yaml_lines(self, cell):
        lines = []
        if cell:
            lines += cell.get_yaml_lines()
        if self._energy is not None:
            lines.append("energy: %20.10f" % self._energy)

        if self._forces is not None:
            lines.append("forces:")
            for i, v in enumerate(self._forces):
                lines.append("- [ %15.10f, %15.10f, %15.10f ] # %d" %
                        (v[0], v[1], v[2], i + 1))

        if self._stress is not None:
            lines.append("stress:")
            for x, v in zip(('x', 'y', 'z'), self._stress):
                lines.append("- [ %15.10f, %15.10f, %15.10f ] # %sx %sy %sz"
                             % (v[0], v[1], v[2], x, x, x))
        return lines

class OneShotCalculation(TaskElement, OneShotCalculationYaml):
    def __init__(self,
                 directory=None,
                 name=None,
                 traverse=False):

        TaskElement.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._tasks = [] # Means singular task
        self._traverse = traverse

        self._cell = None

        self._energy = None
        self._forces = None
        self._stress = None

    def next(self):
        self._collect()
        self._write_yaml()
        raise StopIteration

    def get_cell(self):
        return self._cell

    def get_stress(self):
        return self._stress

    def get_forces(self):
        return self._forces

    def get_energy(self):
        return self._energy

    def set_status(self):
        if self._traverse is False:
            self._status = self._job.get_status()
        else:
            self._status = "done"

    def begin(self):
        if self._traverse is False:
            if self._job is None:
                print "set_job has to be executed."
                raise
            self._prepare() # When traverse != False, files are not created.

        self._status = "begin"

    def done(self):
        return (self._status == "terminate" or
                self._status == "done")

class ElectronicStructureBase(OneShotCalculation):
    def __init__(self,
                 directory="electronic_structure",
                 name=None,
                 traverse=False):

        OneShotCalculation.__init__(self,
                                    directory=directory,
                                    name=name,
                                    traverse=traverse)

        self._task_type = "electronic_structure"
        self._properties = {}

    def get_properties(self):
        return self._properties

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)

        if 'energies' in self._properties:
            self._energy = self._properties['energies'][-1]
        if 'forces' in self._properties:
            self._forces = self._properties['forces'][-1]
        if 'stress' in self._properties:
            self._stress = self._properties['stress'][-1]

        lines += self._get_oneshot_yaml_lines(self._cell)

        if ('eigenvalues' in self._properties and
            'occupancies' in self._properties):
            for i, (e_spin, o_spin) in enumerate(zip(
                    self._properties['eigenvalues'],
                    self._properties['occupancies'])):

                if e_spin is None or o_spin is None:
                    break

                if (len(self._properties['eigenvalues']) == 2 and
                    len(self._properties['occupancies']) == 2):
                    lines.append("eigenvalues_spin%d:" % (i + 1))
                else:
                    lines.append("eigenvalues:")
                for j, (eigs, occs) in enumerate(zip(e_spin, o_spin)):
                    lines.append("- # %d" % (j + 1))
                    for eig, occ in zip(eigs, occs):
                        lines.append("  - [ %15.10f, %15.10f ]" % (eig, occ))
            
        return lines

        
class StructureOptimizationElementBase(OneShotCalculation):
    def __init__(self,
                 directory="structopt_element",
                 name=None,
                 lattice_tolerance=None,
                 force_tolerance=None,
                 pressure_target=None,
                 stress_tolerance=None,
                 max_increase=None,
                 traverse=False):

        OneShotCalculation.__init__(self,
                                    directory=directory,
                                    name=name,
                                    traverse=traverse)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._task_type = "structopt_element"
        self._tasks = [] # Means singular task

        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._traverse = traverse

        self._current_cell = None

    def get_current_cell(self): # cell under structure optimization
        return self._current_cell

    def done(self):
        return (self._status == "terminate" or
                self._status == "done" or
                self._status == "next")

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        lines += self._get_oneshot_yaml_lines(self._current_cell)
        return lines

class ElasticConstantsElementBase(OneShotCalculation):
    def __init__(self,
                 directory="elastic_constants_element",
                 name=None,
                 traverse=False):

        OneShotCalculation.__init__(self,
                                    directory=directory,
                                    name=name,
                                    traverse=traverse)

        self._task_type = "elastic_constants_element"
        self._elastic_constants = None

    def get_elastic_constants(self):
        return self._elastic_constants

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        w.write("status: %s\n" % self._status)
        if self._elastic_constants is not None:
            if self._cell:
                for line in self._cell.get_yaml_lines():
                    w.write(line + "\n")

            w.write("elastic_constants:\n")
            for v in self._elastic_constants:
                w.write("- [ %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f ]\n" % tuple(v))
        w.close()

        

