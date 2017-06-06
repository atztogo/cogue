#include <Python.h>
#include <numpy/arrayobject.h>
#include "xtalcomp_wrapper.h"

static PyObject * py_xtalcomp(PyObject *self, PyObject *args);

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

static PyMethodDef _xtalcomp_methods[] = {
  {"compare", py_xtalcomp, METH_VARARGS, "Dynamical matrix"},
  {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static int _xtalcomp_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _xtalcomp_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_xtalcomp",
  NULL,
  sizeof(struct module_state),
  _xtalcomp_methods,
  NULL,
  _xtalcomp_traverse,
  _xtalcomp_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__xtalcomp(void)

#else
#define INITERROR return

  void
  init_xtalcomp(void)
#endif
{
  struct module_state *st;
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_xtalcomp", _xtalcomp_methods);
#endif

  if (module == NULL)
    INITERROR;

  st = GETSTATE(module);

  st->error = PyErr_NewException("_xtalcomp.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

static PyObject * py_xtalcomp(PyObject *self, PyObject *args)
{
  double tolerance, angle_tolerance;
  PyArrayObject *pylattice1, *pylattice2;
  PyArrayObject *pypositions1, *pypositions2;
  PyArrayObject *pyatom_types1, *pyatom_types2;
  if (!PyArg_ParseTuple(args, "OOOOOOdd",
			&pylattice1,
			&pyatom_types1,
			&pypositions1,
			&pylattice2,
			&pyatom_types2,
			&pypositions2,
			&tolerance,
			&angle_tolerance)) {
    return NULL;
  }

  double (*lattice1)[3] = (double(*)[3])pylattice1->data;
  const int* types1 = (int*)pyatom_types1->data;
  double (*positions1)[3] = (double(*)[3])pypositions1->data;

  double (*lattice2)[3] = (double(*)[3])pylattice2->data;
  const int* types2 = (int*)pyatom_types2->data;
  double (*positions2)[3] = (double(*)[3])pypositions2->data;

  const int num_atom = pyatom_types1->dimensions[0];

  int result = xtalcomp(num_atom,
			lattice1,
			types1,
			positions1,
			lattice2,
			types2,
			positions2,
			tolerance,
			angle_tolerance);
	       
  return PyLong_FromLong((long) result);
}
