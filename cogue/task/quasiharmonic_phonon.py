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

    Four stages:
    0. Structure optimization of input cell
    1. (optional) Equation of state calculation, i.e.,
       total energy calculations with  -2, -1, 0, +1, +2, +3, +4% volumes
    2. (optional) Mode Gruneisen parameter calculations with 0, +1% volumes
    3. Phonon calculations at specified strains or estimated volumes
    
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
        if self._strains is None:
            self._estimate_strain = True
        else:
            self._estimate_strain = False
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

        self._imaginary_ratio = None
        self._grid_spacing = 200

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
            if (task.get_status() == "terminate" or
                task.get_status() == "max_iteration"):
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
            if self._estimate_strain:
                self._all_tasks = [None]
                self._prepare_electric_eos()
            else:
                self._all_tasks = [None, None]
                self._prepare_phonons()
        else:
            self._set_stage0()

    def done(self):
        return (self._status == "done" or
                self._status == "terminate" or
                self._status == "imaginary_modes" or
                self._status == "phonon_for_gruneisen_failed" or
                self._status == "strain_estimation_difficulty" or
                self._status == "next")

    def next(self):
        if self._stage == 0: # Equilibrium
            if self._status == "next":
                if self._estimate_strain:
                    self._prepare_electric_eos()
                    return self._tasks
                else:
                    self._prepare_phonons()
                    return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._set_stage0()
                return self._tasks
        elif self._stage == 1: # EoS
            if self._status == "next":
                if self._estimate_strain:
                    self._prepare_mode_gruneisen()
                    return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                self._prepare_electric_eos()
                return self._tasks
        elif self._stage == 2: # Mode Gruneisen
            if self._status == "next":
                self._prepare_phonons()
                if self._tasks:
                    return self._tasks
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)
                cells = [task.get_cell()
                         for task in self._all_tasks[1].get_all_tasks()[2:5]]
                tasks = self._get_phonon_tasks(cells, is_cell_relaxed=True)
                self._tasks = []
                for i in terminated:
                    self._tasks.append(tasks[i])
                    self._all_tasks[i + 2] = tasks[i]
                self._status = "mode_gruneisen"
                return self._tasks
        else: # QHA
            if self._status == "next":
                if self._calculate_quasiharmonic_phonon():
                    self._status = "done"
                else:
                    self._status = "terminate"
            elif (self._status == "terminate" and self._traverse == "restart"):
                self._traverse = False
                terminated = []
                for i, task in enumerate(self._tasks):
                    if task.get_status() == "terminate":
                        terminated.append(i)

                if len(terminated) > 0:
                    self._tasks = []
                    cells = get_strained_cells(self.get_cell(), self._strains)
                    tasks = self._get_phonon_tasks(cells)
                    for i in terminated:
                        self._tasks.append(tasks[i])
                        if self._estimate_strain:
                            self._all_tasks[i + 5] = tasks[i]
                        else:
                            self._all_tasks[i + 2] = tasks[i]
                    self._status = "phonons"
                    return self._tasks
                
        self._tasks = []
        self._write_yaml()
        raise StopIteration

    def _calculate_quasiharmonic_phonon(self):
        energies = []
        volumes = []
        phonons = []

        if self._estimate_strain:
            phonon_tasks = self._all_tasks[5:]
        else:
            phonon_tasks = self._all_tasks[2:]

        for task in phonon_tasks:
            energies.append(task.get_energy())
            volumes.append(task.get_cell().get_volume())
            phonons.append(task.get_phonon())

        with open("e-v.dat", 'w') as w:
            w.write("#   cell volume        energy of cell other than phonon\n")
            for e, v in zip(energies, volumes):
                w.write("%20.13f %20.13f\n" % (v, e))

        if self._sampling_mesh is not None:
            thermal_properties = []
            for phonon in phonons:
                if not phonon.set_mesh(self._sampling_mesh,
                                       is_gamma_center=self._is_gamma_center):
                    self._log += ("Phonon calculation failed in QHA "
                                  "calculation.\n")
                    self._status = "phonon_for_qha_failed"
                    return False


            for i, phonon in enumerate(phonons):
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

            if qha is None:
                return False
            else:
                qha.write_helmholtz_volume()
                qha.write_volume_temperature()
                qha.write_thermal_expansion()
                qha.write_volume_expansion()
                qha.write_gibbs_temperature()
                qha.write_bulk_modulus_temperature()
                qha.write_heat_capacity_P_numerical()
                qha.write_heat_capacity_P_polyfit()
                qha.write_gruneisen_temperature()
                return True

    def _get_quasiharmonic_phonon(self,
                                  energies,
                                  volumes,
                                  thermal_properties,
                                  t_max):
        T = []
        F = []
        S = []
        Cv = []
        V = []
        U = []

        for tp, v, u in zip(thermal_properties,
                            volumes,
                            energies):
            (temperatures,
             free_energies,
             entropies,
             heat_capacities) = tp

            if (np.isnan(free_energies).any() or
                np.isnan(entropies).any() or
                np.isnan(heat_capacities).any()):
                self._log += "nan is found in thermal property.\n"
                continue

            T.append(temperatures)
            F.append(free_energies)
            S.append(entropies)
            Cv.append(heat_capacities)
            V.append(v)
            U.append(u)

        self._log += "Number of QHA volume points is %d.\n" % len(U)

        if len(U) > 4:
            qha = PhonopyQHA(V,
                             U,
                             temperatures=T[0],
                             free_energy=np.transpose(F),
                             cv=np.transpose(Cv),
                             entropy=np.transpose(S),
                             t_max=t_max,
                             verbose=False)
            return qha
        else:
            return None

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
                 for task in self._all_tasks[1].get_all_tasks()[2:5]]
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

        if self._estimate_strain:
            self._strains = self._get_estimated_strains()

        cell = self.get_cell()
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

    def _check_imaginary(self, phonon, lattice):
        _, weights, freqs, _ = phonon.get_mesh()
        ratio = (float(np.extract(freqs[:, 0] < 0, weights).sum()) /
                 np.prod(self._sampling_mesh))
        return ratio

    def _get_estimated_strains(self):
        t_max = 1500
        t_step = 10
        t_min = 0

        cell = self.get_cell()
        lattice = cell.get_lattice()

        if self._sampling_mesh is None:
            self._sampling_mesh = klength2mesh(self._grid_spacing, lattice)
            self._is_gamma_center = True

        phonons = [task.get_phonon() for task in self._all_tasks[2:5]]
        if None in phonons:
            self._log += ("Phonon calculation failed in mode Gruneisen "
                          "parameter calculation.\n")
            self._status = "phonon_for_gruneisen_failed"
            return []

        if phonons[0].set_mesh(self._sampling_mesh, self._is_gamma_center):
            self._imaginary_ratio = self._check_imaginary(phonons[0], lattice)
        else:
            self._log += ("Phonon calculation failed in mode Gruneisen "
                          "parameter calculation.\n")
            self._status = "phonon_for_gruneisen_failed"
            return []

        if self._imaginary_ratio > 0.01:
            self._log += ("Imaginary modes are found in one point phonon "
                          "calculation.\n")
            self._status = "imaginary_modes"
            return []

        gruneisen = PhonopyGruneisen(phonons[1], phonons[0], phonons[2])

        if not gruneisen.set_mesh(self._sampling_mesh,
                                  is_gamma_center=self._is_gamma_center):
            self._log += ("Phonon calculation failed in mode Gruneisen "
                          "parameter calculation.\n")
            self._status = "phonon_for_gruneisen_failed"
            return []
            
        vol = np.linalg.det(lattice)
        volumes = [vol * (1 + strain) for strain in _eos_strains]
        eos = self._all_tasks[1].get_equation_of_state()
        energies = [eos(v) for v in volumes]

        with open("estimated_e-v.dat", 'w') as w:
            w.write("#   cell volume        energy of cell "
                    "other than phonon\n")
            for e, v in zip(energies, volumes):
                w.write("%20.13f %20.13f\n" % (v, e))

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

        if qha is None:
            self._log += ("Approximated QHA from mode Grunsein parameter "
                          "failed.\n")
            self._status = "strain_estimation_difficulty"
            return []

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
        if self._imaginary_ratio is not None:
            lines.append("imaginary_ratio: %f" % self._imaginary_ratio)
        if self._is_cell_relaxed:
            cell = self._cell
        else:
            cell = self._all_tasks[0].get_cell()
        lines += self._get_phonon_yaml_lines(cell)

        return lines
