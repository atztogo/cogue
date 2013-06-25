.. _VASP_task:

VASP tasks
===========

.. _VASP_electronic_structure_task:

Electronic structure
--------------------------

.. code-block:: python

   def electronic_structure(directory="electronic_structure",
                            name=None,
                            job=None,
                            traverse=False,
                            cell=None,
                            pseudo_potential_map=None,
                            k_mesh=None,
                            k_shift=None,
                            k_gamma=None,
                            k_length=None,
                            incar=None):

.. _VASP_structure_optimization_task:

Structure optimization
----------------------------

.. code-block:: python

   def structure_optimization(directory="structure_optimization",
                              name=None,
                              job=None,
                              lattice_tolerance=0.1,
                              force_tolerance=1e-3,
                              pressure_target=0,
                              stress_tolerance=0.1, # GPa=kbar / 10
                              max_increase=1.5,
                              max_iteration=5,
                              min_iteration=1,
                              find_symmetry=True,
                              impose_symmetry=False,
                              symmetry_tolerance=0.1,
                              traverse=False,
                              cell=None,
                              pseudo_potential_map=None,
                              k_mesh=None,
                              k_shift=None,
                              k_gamma=None,
                              k_length=None,
                              incar=None):

.. _VASP_bulk_modulus_task:

Bulk modulus
------------------

.. code-block:: python

   def bulk_modulus(directory="bulk_modulus",
                    name=None,
                    job=None,
                    lattice_tolerance=0.1,
                    force_tolerance=1e-3,
                    pressure_target=0,
                    stress_tolerance=10,
                    max_increase=1.5,
                    max_iteration=4,
                    min_iteration=1,
                    traverse=False,
                    is_cell_relaxed=False,
                    cell=None,
                    pseudo_potential_map=None,
                    k_mesh=None,
                    k_shift=None,
                    k_gamma=None,
                    k_length=None,
                    incar=None):


Elastic constants
-----------------------

.. _VASP_elastic_constants_task:

.. code-block:: python

   def elastic_constants(directory="elastic_constants",
                         name=None,
                         job=None,
                         lattice_tolerance=0.1,
                         force_tolerance=1e-3,
                         pressure_target=0,
                         stress_tolerance=10,
                         max_increase=1.5,
                         max_iteration=4,
                         min_iteration=1,
                         traverse=False,
                         is_cell_relaxed=False,
                         cell=None,
                         pseudo_potential_map=None,
                         k_mesh=None,
                         k_shift=None,
                         k_gamma=None,
                         k_length=None,
                         incar=None):


.. _VASP_phonon_task:

Phonon
------------

.. code-block:: python

   def phonon(directory="phonon",
              name=None,
              job=None,
              supercell_matrix=np.eye(3, dtype=int),
              primitive_matrix=np.eye(3, dtype=int),
              distance=0.01,
              lattice_tolerance=0.1,
              force_tolerance=1e-3,
              pressure_target=0,
              stress_tolerance=10,
              max_increase=1.5,
              max_iteration=4,
              min_iteration=1,
              traverse=False,
              is_cell_relaxed=False,
              cell=None,
              pseudo_potential_map=None,
              k_mesh=None,
              k_shift=None,
              k_gamma=None,
              k_length=None,
              incar=None):


.. _VASP_mode_gruneisen_task:

Mode-Grüneisen parameter
------------------------------

.. code-block:: python

   def mode_gruneisen(directory="mode_gruneisen",
                      name=None,
                      job=None,
                      delta_strain=0.001,
                      strain=None,
                      bias=None,
                      supercell_matrix=np.eye(3, dtype=int),
                      primitive_matrix=np.eye(3, dtype=int),
                      distance=0.01,
                      lattice_tolerance=0.1,
                      force_tolerance=1e-3,
                      pressure_target=0,
                      stress_tolerance=10,
                      max_increase=1.5,
                      max_iteration=3,
                      min_iteration=1,
                      traverse=False,
                      is_cell_relaxed=False,
                      cell=None,
                      pseudo_potential_map=None,
                      k_mesh=None,
                      k_shift=None,
                      k_gamma=None,
                      k_length=None,
                      incar=None):

.. _VASP_QHA_thermal_expansion_task:

Thermal expansion
-----------------------

.. code-block:: python

   def quasiharmonic_phonon(directory="quasiharmonic_phonon",
                            name=None,
                            job=None,
                            strains=[-0.04, -0.02, 0.02, 0.04, 0.06, 0.08],
                            sampling_mesh=None,
                            t_step=None,
                            t_max=None,
                            t_min=None,
                            supercell_matrix=np.eye(3, dtype=int),
                            primitive_matrix=np.eye(3, dtype=int),
                            distance=0.01,
                            lattice_tolerance=0.1,
                            force_tolerance=1e-3,
                            pressure_target=0,
                            stress_tolerance=10,
                            max_increase=1.5,
                            max_iteration=3,
                            min_iteration=1,
                            traverse=False,
                            is_cell_relaxed=False,
                            cell=None,
                            pseudo_potential_map=None,
                            k_mesh=None,
                            k_shift=None,
                            k_gamma=None,
                            k_length=None,
                            incar=None):

.. _VASP_phonon_relax_task:

Phonon relax
------------------

.. code-block:: python

   def phonon_relax(directory="phonon_relax",
                    name=None,
                    job=None,
                    distance=0.01,
                    lattice_tolerance=0.1,
                    force_tolerance=1e-3,
                    pressure_target=0,
                    stress_tolerance=10,
                    max_increase=1.5,
                    max_iteration=4,
                    min_iteration=1,
                    symmetry_tolerance=0.1,
                    restrict_offspring=False,
                    max_offspring=None,
                    cutoff_eigenvalue=-0.02,
                    max_displacement=None,
                    num_sampling_points=60,
                    traverse=False,
                    cell=None,
                    pseudo_potential_map=None,
                    k_mesh=None,
                    k_shift=None,
                    k_gamma=None,
                    k_length=None,
                    incar=None):

