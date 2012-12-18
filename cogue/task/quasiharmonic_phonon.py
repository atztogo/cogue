from cogue.task import TaskElement
import numpy as np
from phonopy import PhonopyQHA

class QuasiHarmonicPhononBase(TaskElement):
    """QuasiHarmonicPhonon class

    Three stages:
    1. structure optimization of input cell
    2. create cells with various volumes and optimize them
    3. calculate phonons for three cells
    4. calculate quasi-harmonic phonons
    
    """
    
    def __init__(self,
                 directory=None,
                 name=None,
                 strains=None,
                 sampling_mesh=None,
                 t_step=None,
                 t_max=None,
                 t_min=None,
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
        self._task_type = "quasiharmonic_phonon"

        self._lattices = [np.eye(3)]
        for strain in strains:
            if isinstance(strain, int) or isinstance(strain, float):
                self._lattices.append((1 + strain) ** (1.0 / 3) * np.eye(3))
            else:
                self._lattices.append(np.eye(3) + np.array(strain))
        self._sampling_mesh = sampling_mesh
        if t_step is None:
            self._t_step = 10
        else:
            self._t_step = t_step
        if t_max is None:
            self._t_max = 1500
        else:
            self._t_max = t_max
        if t_min is None:
            self._t_min = 0
        else:
            self._t_min = t_min
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
        self._qh_tasks = None

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
            self._qh_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._qh_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._qh_tasks[0]]

    def end(self):
        self._write_yaml()

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):    
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._qh_tasks[0].get_cell())
            else:
                raise StopIteration
        else:
            if self._status == "next":
                self._calculate_quasiharmonic_phonon()
                self._status = "done"
            raise StopIteration

        return self._tasks

    def _calculate_quasiharmonic_phonon(self):
        energies = []
        volumes = []
        phonons = []
        T = []
        F = []
        S = []
        Cv = []
        for i, task in enumerate(self._tasks):
            energies.append(task.get_energy())
            volumes.append(task.get_equilibrium_cell().get_volume())
            if self._sampling_mesh is not None:
                phonon = task.get_phonon()
                phonon.set_mesh(self._sampling_mesh)
                phonon.set_thermal_properties(t_step=self._t_step,
                                              t_max=self._t_max + self._t_step * 2.5,
                                              t_min=self._t_min)
                (temperatures,
                 free_energies,
                 entropies,
                 heat_capacities) = phonon.get_thermal_properties()
                T.append(temperatures)
                F.append(free_energies)
                S.append(entropies)
                Cv.append(heat_capacities)
                phonon.write_yaml_thermal_properties(
                    "thermal_properties-%02d.yaml" % i)

        if self._sampling_mesh:
            qha = PhonopyQHA(volumes,
                             energies,
                             temperatures=T[0],
                             free_energy=np.transpose(F),
                             cv=np.transpose(Cv),
                             entropy=np.transpose(S),
                             t_max=self._t_max,
                             verbose=False)

            qha.write_helmholtz_volume()
            qha.write_volume_temperature()
            qha.write_thermal_expansion()
            qha.write_volume_expansion()
            qha.write_gibbs_temperature()
            qha.write_bulk_modulus_temperature()
            qha.write_heat_capacity_P_numerical()
            qha.write_heat_capacity_P_polyfit()
            qha.write_gruneisen_temperature()
                
        w = open("e-v.dat", 'w')
        w.write("#   cell volume        energy of cell other than phonon\n")
        for e, v in zip(energies, volumes):
            w.write("%20.13f %20.13f\n" % (v, e))
            
        self._quasiharmonic_phonon = None
        
    def _prepare_next(self, cell):
        self._stage = 1
        self._status = "phonons"
        self._qh_tasks += self._get_phonon_tasks(cell)
        self._tasks = self._qh_tasks[1:]

    def _write_yaml(self):
        w = open("%s.yaml" % self._directory, 'w')
        if self._qh_tasks[0]:
            w.write("lattice_tolerance: %f\n" % self._lattice_tolerance)
            w.write("pressure_target: %f\n" % self._pressure_target)
            w.write("stress_tolerance: %f\n" % self._stress_tolerance)
            w.write("force_tolerance: %f\n" % self._force_tolerance)
            w.write("max_increase: %f\n" % self._max_increase)
            w.write("max_iteration: %d\n" % self._max_iteration)
            w.write("min_iteration: %d\n" % self._min_iteration)
            w.write("iteration: %d\n" % self._qh_tasks[0].get_stage())
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
            cell = self._qh_tasks[0].get_cell()

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
