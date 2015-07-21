#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

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

struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
  struct module_state *st = GETSTATE(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _spglib_methods[] = {
  {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
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

#if PY_MAJOR_VERSION >= 3

static int _spglib_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _spglib_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_spglib",
  NULL,
  sizeof(struct module_state),
  _spglib_methods,
  NULL,
  _spglib_traverse,
  _spglib_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__spglib(void)

#else
#define INITERROR return

  void
  init_spglib(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_spglib", _spglib_methods);
#endif

  if (module == NULL)
    INITERROR;
  struct module_state *st = GETSTATE(module);

  st->error = PyErr_NewException("_spglib.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

static PyObject * get_dataset(PyObject *self, PyObject *args)
{
  int i, j, k, n;
  double symprec;
  SpglibDataset *dataset;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  PyObject* array, *vec, *mat, *rot, *trans, *wyckoffs, *equiv_atoms;
  PyObject *std_lattice, *std_types, *std_positions;
  
  if (!PyArg_ParseTuple(args, "OOOd",
			&lattice,
			&position,
			&atom_type,
			&symprec)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  const int num_atom = PyArray_DIMS(position)[0];
  const int* typat = (int*)PyArray_DATA(atom_type);

  dataset = spg_get_dataset(lat, pos, typat, num_atom, symprec);
  
  array = PyList_New(13);
  n = 0;

  /* Space group number, international symbol, hall symbol */
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->spacegroup_number));
  n++;
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->hall_number));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->international_symbol));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->hall_symbol));
  n++;

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }
  PyList_SetItem(array, n, mat);
  n++;

  /* Origin shift */
  vec = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(vec, i, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_SetItem(array, n, vec);
  n++;

  /* Rotation matrices */
  rot = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(3);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(3);
      for (k = 0; k < 3; k++) {
	PyList_SetItem(vec, k, PyLong_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_SetItem(mat, j, vec);
    }
    PyList_SetItem(rot, i, mat);
  }
  PyList_SetItem(array, n, rot);
  n++;

  /* Translation vectors */
  trans = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_SetItem(trans, i, vec);
  }
  PyList_SetItem(array, n, trans);
  n++;

  /* Wyckoff letters, Equivalent atoms */
  wyckoffs = PyList_New(dataset->n_atoms);
  equiv_atoms = PyList_New(dataset->n_atoms);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_SetItem(wyckoffs, i, PyLong_FromLong((long) dataset->wyckoffs[i]));
    PyList_SetItem(equiv_atoms, i, PyLong_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_SetItem(array, n, wyckoffs);
  n++;
  PyList_SetItem(array, n, equiv_atoms);
  n++;

  std_lattice = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->std_lattice[i][j]));
    }
    PyList_SetItem(std_lattice, i, vec);
  }
  PyList_SetItem(array, n, std_lattice);
  n++;

  std_types = PyList_New(dataset->n_std_atoms);
  std_positions = PyList_New(dataset->n_std_atoms);
  for (i = 0; i < dataset->n_std_atoms; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->std_positions[i][j]));
    }
    PyList_SetItem(std_types, i, PyLong_FromLong((long) dataset->std_types[i]));
    PyList_SetItem(std_positions, i, vec);
  }
  PyList_SetItem(array, n, std_types);
  n++;
  PyList_SetItem(array, n, std_positions);
  n++;

  spg_free_dataset(dataset);

  return array;
}

static PyObject * get_crystallographic_cell(PyObject *self, PyObject *args)
{
  int num_atom, num_atom_std;
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

  p_lattice = (double(*))PyArray_DATA(lattice_vectors);
  p_positions = (double(*))PyArray_DATA(atomic_positions);
  num_atom = PyArray_DIMS(atom_types)[0];
  positions = (double(*)[3]) malloc(sizeof(double[3]) * num_atom * 4);
  types_int = (int*)PyArray_DATA(atom_types);
  types = (int*) malloc(sizeof(int) * num_atom * 4);

  set_spglib_cell(lattice, positions, types, num_atom,
		  p_lattice, p_positions, types_int);
  num_atom_std = spg_refine_cell(lattice, positions, types, num_atom, symprec);
  array = set_cell(num_atom_std, lattice, positions, types);

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

  p_lattice = (double(*))PyArray_DATA(lattice_vectors);
  p_positions = (double(*))PyArray_DATA(atomic_positions);
  types_int = (int*)PyArray_DATA(atom_types);

  num_atom = PyArray_DIMS(atom_types)[0];
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

