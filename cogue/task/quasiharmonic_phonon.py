from cogue.task import TaskElement
from cogue.task.phonon import PhononYaml
from cogue.crystal.cell import get_strained_cells
from cogue.crystal.supercell import estimate_supercell_matrix
from cogue.crystal.utility import klength2mesh
import numpy as np
from phonopy import PhonopyQHA
from phonopy import PhonopyGruneisen

_eos_strains = [-0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04]

class QuasiHarmonicPhononBase(TaskElement, PhononYaml):
    """QuasiHarmonicPhonon class

    Three stages:
    0. Structure optimization of input cell
    1. (optional) Equation of state calculation, i.e.,
       total energy calculations with  -2, -1, 0, +1, +2, +3, +4% volumes
    2. (optional) Mode Gruneisen parameter calculations with 0, +1% volumes
    3. Phonon calculations at specified strains or estimated volumes
    4. Calculate thermal properties
    
    """
    
    def __init__(self,
                 directory=None,
                 name=None,
                 strains=None,
                 sampling_mesh=None,
                 is_gamma_center=False,
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

        self._strains = strains
        self._sampling_mesh = sampling_mesh
        self._is_gamma_center = is_gamma_center
        if t_step is None:
            self._t_step = 2
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

        self._write_yaml()

    def begin(self):
        if not self._job:
            print "set_job has to be executed."
            raise

        if self._is_cell_relaxed:
            if self._strains:
                self._all_tasks = [None, None]
                self._prepare_phonons()
            else:
                self._all_tasks = [None]
                self._prepare_electric_eos()
        else:
            self._set_stage0()

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "next")

    def next(self):
        if self._stage == 0:
            if self._status == "next":
                if self._strains is None:
                    self._prepare_electric_eos()
                    return self._tasks
                else:
                    self._prepare_phonons()
                    return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._set_stage0()
                return self._tasks
        elif self._stage == 1:
            if self._status == "next":
                if self._strains is None:
                    self._prepare_mode_gruneisen()
                    return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._prepare_electric_eos()
                return self._tasks
        elif self._stage == 2:
            if self._status == "next":
                self._prepare_phonons()
                return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)
                cells = [task.get_cell()
                         for task in self._all_tasks[1].get_all_tasks()[2:4]]
                tasks = self._get_phonon_tasks(cells, is_cell_relaxed=True)
                self._tasks = []
                for i in terminated:
                    self._tasks.append(tasks[i])
                    self._all_tasks[i + 2] = tasks[i]
                self._status = "mode_gruneisen"
                return self._tasks
        else:
            if self._status == "next":
                self._calculate_quasiharmonic_phonon()
                self._status = "done"
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)

                if self._strains is None:
                    cells = self._estimate_cells()
                else:
                    cells = get_strained_cells(self.get_cell(), self._strains)
                tasks = self._get_phonon_tasks(cells)
                self._tasks = []
                for i in terminated:
                    self._tasks.append(tasks[i])
                    if self._strains is None:
                        self._all_tasks[i + 4] = tasks[i]
                    else:
                        self._all_tasks[i + 2] = tasks[i]
                self._status = "phonons"
                return self._tasks
                
        self._write_yaml()
        raise StopIteration

    def _calculate_quasiharmonic_phonon(self):
        energies = []
        volumes = []
        phonons = []

        if self._strains is None:
            phonon_tasks = self._all_tasks[4:]
            # energies += [task.get_energy()
            #              for task in self._all_tasks[1].get_all_tasks()[2:4]]
            # volumes += [task.get_cell().get_volume()
            #             for task in self._all_tasks[2:4]]
            # phonons += [task.get_phonon() for task in self._all_tasks[2:4]]
        else:
            phonon_tasks = self._all_tasks[1:]

        for task in phonon_tasks:
            energies.append(task.get_energy())
            volumes.append(task.get_cell().get_volume())
            phonons.append(task.get_phonon())

        if self._sampling_mesh is not None:
            thermal_properties = []
            for i, phonon in enumerate(phonons):
                phonon.set_mesh(self._sampling_mesh,
                                is_gamma_center=self._is_gamma_center)
                phonon.set_thermal_properties(
                    t_step=self._t_step,
                    t_max=self._t_max + self._t_step * 3.5,
                    t_min=self._t_min)
                thermal_properties.append(phonon.get_thermal_properties())
                phonon.write_yaml_thermal_properties(
                    filename="thermal_properties-%02d.yaml" % i)

            qha = self._get_quasiharmonic_phonon(energies,
                                                 volumes,
                                                 thermal_properties,
                                                 self._t_max)

            qha.write_helmholtz_volume()
            qha.write_volume_temperature()
            qha.write_thermal_expansion()
            qha.write_volume_expansion()
            qha.write_gibbs_temperature()
            qha.write_bulk_modulus_temperature()
            qha.write_heat_capacity_P_numerical()
            qha.write_heat_capacity_P_polyfit()
            qha.write_gruneisen_temperature()

        with open("e-v.dat", 'w') as w:
            w.write("#   cell volume        energy of cell other than phonon\n")
            for e, v in zip(energies, volumes):
                w.write("%20.13f %20.13f\n" % (v, e))
            
        self._quasiharmonic_phonon = None

    def _get_quasiharmonic_phonon(self,
                                  energies,
                                  volumes,
                                  thermal_properties,
                                  t_max):
        T = []
        F = []
        S = []
        Cv = []

        for tp in thermal_properties:
            (temperatures,
             free_energies,
             entropies,
             heat_capacities) = tp
            T.append(temperatures)
            F.append(free_energies)
            S.append(entropies)
            Cv.append(heat_capacities)

        qha = PhonopyQHA(volumes,
                         energies,
                         temperatures=T[0],
                         free_energy=np.transpose(F),
                         cv=np.transpose(Cv),
                         entropy=np.transpose(S),
                         t_max=t_max,
                         verbose=False)
        return qha

    def _set_stage0(self):
        self._stage = 0
        self._status = "equilibrium"
        task = self._get_equilibrium_task()
        self._all_tasks = [task]
        self._tasks = [task]

    def _prepare_electric_eos(self):
        self._stage = 1
        self._status = "electric_eos"

        cell = self.get_cell()
        task = self._get_bulk_modulus_task(cell, _eos_strains)
        self._all_tasks.append(task)
        self._tasks = [task]

    def _prepare_mode_gruneisen(self):
        self._stage = 2
        self._status = "mode_gruneisen"

        if self._supercell_matrix is None:
            cell = self.get_cell()
            self._supercell_matrix = estimate_supercell_matrix(
                cell,
                max_num_atoms=self._max_num_atoms)

        cells = [task.get_cell()
                 for task in self._all_tasks[1].get_all_tasks()[2:4]]
        tasks = self._get_phonon_tasks(cells,
                                       is_cell_relaxed=True,
                                       directory="gruneisen")
        self._all_tasks += tasks
        self._tasks = tasks

    def _prepare_phonons(self):
        self._stage = 3
        self._status = "phonons"

        if self._supercell_matrix is None:
            self._supercell_matrix = estimate_supercell_matrix(
                cell,
                max_num_atoms=self._max_num_atoms)

        cell = self.get_cell()
        if self._strains is None:
            cells = self._estimate_cells()
            tasks = self._get_phonon_tasks(cells)
        else:
            cells = get_strained_cells(cell, self._strains)
            tasks = self._get_phonon_tasks(cells)

        self._all_tasks += tasks
        self._tasks = tasks

    def _get_phonon_tasks(self,
                          cells,
                          is_cell_relaxed=False,
                          directory="phonon"):
        phonons = []
        for i, cell in enumerate(cells):
            phonons.append(
                self._get_phonon_task(cell,
                                      "%s-%02d" % (directory, i),
                                      is_cell_relaxed=is_cell_relaxed))

        return phonons

    def _estimate_cells(self):
        strains = self._get_estimated_strains()
        cell = self.get_cell()
        return get_strained_cells(cell, strains)

    def _get_estimated_strains(self, distance=200):
        t_max = 1500
        t_step = 10
        t_min = 0

        phonons = [task.get_phonon() for task in self._all_tasks[-2:]]
        gruneisen = PhonopyGruneisen(phonons[0], phonons[0], phonons[1])
        cell = self.get_cell()
        lattice = cell.get_lattice()

        if self._sampling_mesh is None:
            self._sampling_mesh = klength2mesh(distance, lattice)
            self._is_gamma_center = True

        gruneisen.set_mesh(self._sampling_mesh,
                           is_gamma_center=self._is_gamma_center)
        vol = np.linalg.det(lattice)
        volumes = [vol * (1 + strain) for strain in _eos_strains]
        eos = self._all_tasks[1].get_equation_of_state()
        energies = [eos(v) for v in volumes]
        gruneisen.set_thermal_properties(volumes,
                                         t_step=t_step,
                                         t_max=t_max + t_step * 3.5,
                                         t_min=t_min,
                                         cutoff_frequency=0.1)
        gruneisen.write_yaml_thermal_properties(
            filename="estimated_thermal_props")
        gruneisen_tp = gruneisen.get_thermal_properties()
        thermal_properties = [tp.get_thermal_properties()
                              for tp in gruneisen_tp.get_thermal_properties()]
        qha = self._get_quasiharmonic_phonon(energies,
                                             volumes,
                                             thermal_properties,
                                             t_max)

        equi_volumes = qha.get_volume_temperature()
        #           0K 1000K
        # |---|---|-o-|-o-|---|---|---|---|---|
        # 0   1   2   3   4   5   6   7   8   9
        #
        d_v = (equi_volumes[100] - equi_volumes[0])
        d_left = equi_volumes[0] - d_v * 2.5
        d_right = equi_volumes[100] + d_v * 5.5
        strains = [v / vol - 1 for v in np.linspace(d_left, d_right, 10)]
        return strains

    def get_yaml_lines(self):
        lines = TaskElement.get_yaml_lines(self)
        if self._sampling_mesh is not None:
            lines.append("sampling_mesh: [ %3d, %3d, %3d ]" %
                         tuple(self._sampling_mesh))
            if self._is_gamma_center:
                lines.append("is_gamma_center: True")
            else:
                lines.append("is_gamma_center: False")
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._all_tasks[0].get_cell()
        lines += self._get_phonon_yaml_lines(cell)

        return lines
