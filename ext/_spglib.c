#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_dataset(PyObject *self, PyObject *args);
static PyObject * get_crystallographic_cell(PyObject *self, PyObject *args);
/* static PyObject * find_primitive(PyObject *self, PyObject *args); */

static PyMethodDef functions[] = {
  {"get_dataset", get_dataset, METH_VARARGS,
   "Return crystal symmetry dataset"},
  {"get_crystallographic_cell", get_crystallographic_cell, METH_VARARGS,
   "Return crystallographically proper cell"},
  /* {"find_primitive", find_primitive, METH_VARARGS, "Find primitive cell in the input cell"}, */
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
  
  if (!PyArg_ParseTuple(args, "OOOd",
			 &lattice_vectors,
			 &atomic_positions,
			 &atom_types,
			 &symprec)) {
    return NULL;
  }

  double * p_lattice = (double(*))lattice_vectors->data;
  double * p_positions = (double(*))atomic_positions->data;
  double lattice[3][3];
  num_atom = atomic_positions->dimensions[1];
  double (*positions)[3];
  positions = (double(*) [3]) malloc(sizeof(double[3]) * num_atom);
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice[i][j] = p_lattice[ i * 3 + j ];
    }
  }

  for (i = 0; i < num_atom; i++) {
    for (j = 0; j < 3; j++) {
      positions[i][j] = p_positions[ j * num_atom + i ];
    }
  }

  const long* types_long = (long*)atom_types->data;
  int *types;

  types = (int*) malloc(sizeof(int) * num_atom);
  for (i = 0; i < num_atom; i++) {
    types[i] = (int) types_long[i];
  }

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
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
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
	PyList_SetItem(vec, k, PyInt_FromLong((long) dataset->rotations[i][j][k]));
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
    PyList_SetItem(equiv_atoms, i, PyInt_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_SetItem(array, 8, wyckoffs);
  PyList_SetItem(array, 9, equiv_atoms);

  spg_free_dataset(dataset);

  return array;
}

static PyObject * get_crystallographic_cell(PyObject *self, PyObject *args)
{
  int i, j, num_atom, num_atom_brv;
  double symprec;
  PyArrayObject* lattice_vectors;
  PyArrayObject* atomic_positions;
  PyArrayObject* atom_types;
  PyObject *points, *brv_lattice, *numbers, *vec, *array;
  if (!PyArg_ParseTuple(args,
			"OOOd",
			&lattice_vectors,
			&atomic_positions,
			&atom_types,
			&symprec)) {
    return NULL;
  }


  double * p_lattice = (double(*))lattice_vectors->data;
  double * p_positions = (double(*))atomic_positions->data;
  double lattice[3][3];
  num_atom = atomic_positions->dimensions[1];
  double (*positions)[3];
  positions = (double(*) [3]) malloc(sizeof(double[3]) * num_atom * 4);
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice[i][j] = p_lattice[ i * 3 + j ];
    }
  }

  for (i = 0; i < num_atom; i++) {
    for (j = 0; j < 3; j++) {
      positions[i][j] = p_positions[ j * num_atom + i ];
    }
  }

  const long* types_long = (long*)atom_types->data;
  int *types;
  types = (int*) malloc(sizeof(int) * num_atom * 4);

  for (i = 0; i < num_atom; i++) {
    types[i] = (int) types_long[i];
  }

  num_atom_brv = spg_refine_cell(lattice, positions, types, num_atom, symprec);

  array = PyList_New(3);

  /* Lattice */
  brv_lattice = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(lattice[i][j]));
    }
    PyList_SetItem(brv_lattice, i, vec);
  }
  PyList_SetItem(array, 0, brv_lattice);

  /* Points */
  points = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(num_atom_brv);
    for (j = 0; j < num_atom_brv; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(positions[j][i]));
    }
    PyList_SetItem(points, i, vec);
  }
  PyList_SetItem(array, 1, points);

  /* Numbers */
  numbers = PyList_New(num_atom_brv);
  for (i = 0; i < num_atom_brv; i++) {
    PyList_SetItem(numbers, i, PyInt_FromLong((long) types[i]));
  }
  PyList_SetItem(array, 2, numbers);

  free(types);
  free(positions);

  return array;
}


/* static PyObject * find_primitive(PyObject *self, PyObject *args) */
/* { */
/*   int i; */
/*   double symprec; */
/*   PyArrayObject* lattice; */
/*   PyArrayObject* position; */
/*   PyArrayObject* atom_type; */
/*   if (!PyArg_ParseTuple(args, "OOOd", &lattice, &position, &atom_type, &symprec)) */
/*     return NULL; */

/*   double (*lat)[3] = (double(*)[3])lattice->data; */
/*   double (*pos)[3] = (double(*)[3])position->data; */
/*   int num_atom = position->dimensions[0]; */
/*   long* types_long = (long*)atom_type->data; */

/*   int types[num_atom]; */
/*   for (i = 0; i < num_atom; i++) { */
/*     types[i] = (int)types_long[i]; */
/*   } */
  
/*   int num_atom_prim = spg_find_primitive(lat, */
/* 					  pos, */
/* 					  types, */
/* 					  num_atom, */
/* 					  symprec); */

/*   for (i = 0; i < num_atom_prim; i++) { */
/*     types_long[i] = types[i]; */
/*   } */

/*   return PyInt_FromLong((long) num_atom_prim); */
/* } */

