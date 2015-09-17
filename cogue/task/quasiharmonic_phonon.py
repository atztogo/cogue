from cogue.task import TaskElement
from cogue.task.phonon import PhononYaml
from cogue.crystal.supercell import estimate_supercell_matrix
import numpy as np
from phonopy import PhonopyQHA

class QuasiHarmonicPhononBase(TaskElement, PhononYaml):
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
                 is_cell_relaxed=False,
                 max_num_atoms=None,
                 traverse=False):

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
        self._max_num_atoms = max_num_atoms
        
        self._stage = 0
        self._tasks = None

        self._cell = None
        self._all_tasks = None

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

        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            self._all_tasks = [None]
            self._prepare_next(self._cell)
        else:
            self._status = "equilibrium"
            self._all_tasks = [self._get_equilibrium_task()]
            self._tasks = [self._all_tasks[0]]

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):    
        if self._stage == 0:
            if self._status == "next":
                self._prepare_next(self._all_tasks[0].get_cell())
                return self._tasks
        else:
            if self._status == "next":
                self._calculate_quasiharmonic_phonon()
                self._status = "done"
                
        self._write_yaml()
        raise StopIteration

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
            volumes.append(task.get_cell().get_volume())
            if self._sampling_mesh is not None:
                phonon = task.get_phonon()
                phonon.set_mesh(self._sampling_mesh)
                phonon.set_thermal_properties(
                    t_step=self._t_step,
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

        if self._supercell_matrix is None:
            self._supercell_matrix = estimate_supercell_matrix(
                cell,
                max_num_atoms=self._max_num_atoms)

        self._all_tasks += self._get_phonon_tasks(cell)
        self._tasks = self._all_tasks[1:]

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._all_tasks[0].get_cell()
        lines += self._get_phonon_yaml_lines(cell)

        return lines

