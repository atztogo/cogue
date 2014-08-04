#include <Python.h>
#include <numpy/arrayobject.h>
#include "xtalcomp_wrapper.h"

static PyObject * py_xtalcomp(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"compare", py_xtalcomp, METH_VARARGS, "Dynamical matrix"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_xtalcomp(void)
{
  Py_InitModule3("_xtalcomp", functions, "C-extension for xtalcomp\n\n...\n");
  return;
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
	       
  return PyInt_FromLong((long) result);
}
