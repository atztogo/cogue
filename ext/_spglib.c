#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_dataset(PyObject *self, PyObject *args);
static PyObject * get_crystallographic_cell(PyObject *self, PyObject *args);
static PyObject * get_primitive(PyObject *self, PyObject *args);
static PyObject * get_datasets_of_modulations(PyObject *self, PyObject *args);
static void set_spglib_cell(double lattice[3][3],
			    double positions[][3],
			    int types[],
			    const int num_atom,
			    const double * p_lattice,
			    const double * p_positions,
			    const int * types_int);
static PyObject * set_cell(int num_atom,
			   double lattice[3][3],
			   double positions[][3],
			   int types[]);

static PyMethodDef functions[] = {
  {"get_dataset", get_dataset, METH_VARARGS,
   "Return crystal symmetry dataset"},
  {"get_crystallographic_cell", get_crystallographic_cell, METH_VARARGS,
   "Return crystallographically proper cell"},
  {"get_primitive_cell", get_primitive, METH_VARARGS,
   "Return a primitive cell"},
  {"get_datasets_of_modulations", get_datasets_of_modulations, METH_VARARGS,
   "Return datasets of modulations"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_spglib(void)
{
  Py_InitModule3("_spglib", functions, "C-extension for spglib\n\n...\n");
  return;
}

static PyObject * get_dataset(PyObject *self, PyObject *args)
{
  int i, j, k, num_atom;
  double symprec;
  SpglibDataset *dataset;
  PyArrayObject* lattice_vectors;
  PyArrayObject* atomic_positions;
  PyArrayObject* atom_types;
  PyObject* array, *vec, *mat, *rot, *trans, *wyckoffs, *equiv_atoms;
  
  double *p_lattice;
  double *p_positions;
  double lattice[3][3];
  double (*positions)[3];
  int *types_int;
  int *types;

  if (!PyArg_ParseTuple(args, "OOOd",
			&lattice_vectors,
			&atomic_positions,
			&atom_types,
			&symprec)) {
    return NULL;
  }

  p_lattice = (double(*))lattice_vectors->data;
  p_positions = (double(*))atomic_positions->data;
  num_atom = atom_types->dimensions[0];
  positions = (double(*)[3]) malloc(sizeof(double[3]) * num_atom);
  types_int = (int*)atom_types->data;
  types = (int*) malloc(sizeof(int) * num_atom);
  set_spglib_cell(lattice, positions, types, num_atom,
		  p_lattice, p_positions, types_int);
  dataset = spg_get_dataset(lattice, positions, types, num_atom, symprec);
  free(types);
  free(positions);

  array = PyList_New(10);

  /* Space group number, international symbol, hall symbol */
  PyList_SetItem(array, 0, PyInt_FromLong((long) dataset->spacegroup_number));
  PyList_SetItem(array, 1, PyString_FromString(dataset->international_symbol));
  PyList_SetItem(array, 2, PyInt_FromLong((long) dataset->hall_number));
  PyList_SetItem(array, 3, PyString_FromString(dataset->hall_symbol));

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j,
		     PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }
  PyList_SetItem(array, 4, mat);

  /* Origin shift */
  vec = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(vec, i, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_SetItem(array, 5, vec);

  /* Rotation matrices */
  rot = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(3);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(3);
      for (k = 0; k < 3; k++) {
	PyList_SetItem(vec, k,
		       PyInt_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_SetItem(mat, j, vec);
    }
    PyList_SetItem(rot, i, mat);
  }
  PyList_SetItem(array, 6, rot);

  /* Translation vectors */
  trans = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_SetItem(trans, i, vec);
  }
  PyList_SetItem(array, 7, trans);

  /* Wyckoff letters, Equivalent atoms */
  wyckoffs = PyList_New(dataset->n_atoms);
  equiv_atoms = PyList_New(dataset->n_atoms);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_SetItem(wyckoffs, i, PyInt_FromLong((long) dataset->wyckoffs[i]));
    PyList_SetItem(equiv_atoms, i,
		   PyInt_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_SetItem(array, 8, wyckoffs);
  PyList_SetItem(array, 9, equiv_atoms);

  spg_free_dataset(dataset);

  return array;
}

static PyObject * get_crystallographic_cell(PyObject *self, PyObject *args)
{
  int num_atom, num_atom_brv;
  double symprec;
  PyArrayObject* lattice_vectors;
  PyArrayObject* atomic_positions;
  PyArrayObject* atom_types;
  PyObject *array;

  double *p_lattice;
  double *p_positions;
  double lattice[3][3];
  double (*positions)[3];
  int *types_int;
  int *types;

  if (!PyArg_ParseTuple(args,
			"OOOd",
			&lattice_vectors,
			&atomic_positions,
			&atom_types,
			&symprec)) {
    return NULL;
  }

  p_lattice = (double(*))lattice_vectors->data;
  p_positions = (double(*))atomic_positions->data;
  num_atom = atom_types->dimensions[0];
  positions = (double(*)[3]) malloc(sizeof(double[3]) * num_atom * 4);
  types_int = (int*)atom_types->data;
  types = (int*) malloc(sizeof(int) * num_atom * 4);

  set_spglib_cell(lattice, positions, types, num_atom,
		  p_lattice, p_positions, types_int);
  num_atom_brv = spg_refine_cell(lattice, positions, types, num_atom, symprec);
  array = set_cell(num_atom_brv, lattice, positions, types);

  free(types);
  free(positions);

  return array;
}


static PyObject * get_primitive(PyObject *self, PyObject *args)
{
  int num_atom, num_atom_prim;
  double symprec;
  PyArrayObject* lattice_vectors;
  PyArrayObject* atomic_positions;
  PyArrayObject* atom_types;
  PyObject *array;

  double *p_lattice;
  double *p_positions;
  double lattice[3][3];
  double (*positions)[3];
  int *types_int;
  int *types;

  if (!PyArg_ParseTuple(args,
			"OOOd",
			&lattice_vectors,
			&atomic_positions,
			&atom_types,
			&symprec)) {
    return NULL;
  }

  p_lattice = (double(*))lattice_vectors->data;
  p_positions = (double(*))atomic_positions->data;
  types_int = (int*)atom_types->data;

  num_atom = atom_types->dimensions[0];
  positions = (double(*)[3]) malloc(sizeof(double[3]) * num_atom);
  types = (int*) malloc(sizeof(int) * num_atom);
  set_spglib_cell(lattice, positions, types, num_atom,
		  p_lattice, p_positions, types_int);
  num_atom_prim = spg_find_primitive(lattice, positions, types, num_atom,
				     symprec);
  array = set_cell(num_atom_prim, lattice, positions, types);

  free(types);
  free(positions);

  return array;
}

static PyObject * get_datasets_of_modulations(PyObject *self, PyObject *args)
{
  /* This function aims improving the speed of symmetry searches of */
  /* huge configurations of combinations of degenerated modulations. */
  /* The plan: */
  /* 1. Recieve a set of modulation basis (e.g., {v_a, v_b, v_c}) and those */
  /*    configurations (e.g., {(1, 0, 0), (sqrt(2)/2, 0, sqrt(2)/2), ...} */
  /*    norm = 1). */
  /* 2. Create real value modulations. */
  /* 3. Shift supercell points by the modulations. */
  /* 4. Search symmetry of the modified supercell. */
  /* These operations may be multithreaded. */
}

static PyObject * set_cell(int num_atom,
			   double lattice[3][3],
			   double positions[][3],
			   int types[])
{
  int i, j;
  PyObject *points, *new_lattice, *numbers, *vec, *array;

  array = PyList_New(3);

  /* Lattice */
  new_lattice = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(lattice[i][j]));
    }
    PyList_SetItem(new_lattice, i, vec);
  }
  PyList_SetItem(array, 0, new_lattice);

  /* Points */
  points = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(num_atom);
    for (j = 0; j < num_atom; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(positions[j][i]));
    }
    PyList_SetItem(points, i, vec);
  }
  PyList_SetItem(array, 1, points);

  /* Numbers */
  numbers = PyList_New(num_atom);
  for (i = 0; i < num_atom; i++) {
    PyList_SetItem(numbers, i, PyInt_FromLong((long) types[i]));
  }
  PyList_SetItem(array, 2, numbers);

  return array;
}

static void set_spglib_cell(double lattice[3][3],
			    double positions[][3],
			    int types[],
			    const int num_atom,
			    const double * p_lattice,
			    const double * p_positions,
			    const int * types_int)
{
  int i, j;
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice[i][j] = p_lattice[i * 3 + j];
    }
  }

  for (i = 0; i < num_atom; i++) {
    for (j = 0; j < 3; j++) {
      positions[i][j] = p_positions[j * num_atom + i];
    }
  }

  for (i = 0; i < num_atom; i++) {
    types[i] = (int) types_int[i];
  }
}

