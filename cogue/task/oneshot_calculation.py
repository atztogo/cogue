import os
from cogue.task import TaskElement

class OneShotCalculation(TaskElement):
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

    def get_cell(self):
        return self._cell

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

    def end(self):
        self._collect()
        self._write_yaml()

    def done(self):
        return ("terminate" in self._status or 
                "done" in self._status)



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

    def _write_yaml(self):
        w = open("electronic_structure.yaml", 'w')
        w.write("status: %s\n" % self._status)
        if 'energies' in self._properties:
            w.write("energy: %20.10f\n" % self._properties['energies'][-1])
        if 'forces' in self._properties:
            w.write("forces:\n")
            for i, v in enumerate(self._properties['forces'][-1]):
                w.write("- [ %15.10f, %15.10f, %15.10f ] # %d\n" %
                        (v[0], v[1], v[2], i + 1))
        if 'stress' in self._properties:
            w.write("stress:\n")
            for v in self._properties['stress'][-1]:
                w.write("- [ %15.10f, %15.10f, %15.10f ]\n" %
                        tuple(v))
        if 'eigenvalues' in self._properties:
            if len(self._properties['eigenvalues']) == 2:
                w.write("eigenvalues_spin1:\n")
            else:
                w.write("eigenvalues:\n")
            for i, (eigs, occs) in enumerate(zip(
                    self._properties['eigenvalues'][0],
                    self._properties['occupancies'][0])):
                w.write("- # %d\n" % (i + 1))
                for eig, occ in zip(eigs, occs):
                    w.write("  - [ %15.10f, %15.10f ]\n" % (eig, occ))
            if len(self._properties['eigenvalues']) == 2:
                w.write("eigenvalues_spin2:\n")
                for i, (eigs, occs) in enumerate(zip(
                        self._properties['eigenvalues'][1],
                        self._properties['occupancies'][1])):
                    w.write("- # %d\n" % (i + 1))
                    for eig, occ in zip(eigs, occs):
                        w.write("  - [ %15.10f, %15.10f ]\n" % (eig, occ))
            
        w.close()

        
class StructureOptimizationElementBase(OneShotCalculation):
    def __init__(self,
                 directory="structure_optimization_element",
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
        self._task_type = "structure_optimization_element"
        self._tasks = [] # Means singular task

        self._lattice_tolerance = lattice_tolerance
        self._pressure_target = pressure_target
        self._stress_tolerance = stress_tolerance
        self._force_tolerance = force_tolerance
        self._max_increase = max_increase
        self._traverse = traverse

        self._forces = None
        self._stress = None
        self._energy = None

        self._current_cell = None

    def get_current_cell(self): # cell under structure optimization
        return self._current_cell

    def get_stress(self):
        return self._stress

    def get_forces(self):
        return self._forces

    def get_energy(self):
        return self._energy

    def done(self):
        return ("terminate" in self._status or
                "done" in self._status or
                "next" in self._status)

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        w.write("status: %s\n" % self._status)

        if self._current_cell:
            lattice = self._current_cell.get_lattice().T
            points = self._current_cell.get_points().T
            symbols = self._current_cell.get_symbols()
        
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
        if not self._elastic_constants == None:
            w.write("elastic_constants:\n")
            for v in self._elastic_constants:
                w.write("- [ %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f ]\n" % tuple(v))
        w.close()

        

