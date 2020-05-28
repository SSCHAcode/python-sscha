#include <python2.7/Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "LanczosFunctions.h"

#ifdef _MPI
#include <mpi.h>
#endif

#include <omp.h>

static PyObject *GetV3 (PyObject * self, PyObject * args);
static PyObject *GetV4 (PyObject * self, PyObject * args);
static PyObject *ApplyV3ToDyn(PyObject * self, PyObject * args);
static PyObject *ApplyV3ToVector(PyObject * self, PyObject * args);
static PyObject *ApplyV4ToDyn(PyObject * self, PyObject * args);
static PyObject *ApplyV3_FT(PyObject * self, PyObject * args);
static PyObject *ApplyV4_FT(PyObject * self, PyObject * args);


static PyMethodDef odd_engine[] = {
    {"GetV3", GetV3, METH_VARARGS, "Compute the v3 matrix using the replica method"},
    {"GetV4", GetV4, METH_VARARGS, "Compute the sparse v4 using the replica method"},
    {"ApplyV3ToDyn", ApplyV3ToDyn, METH_VARARGS, "Apply the v3 to a given dynamical matrix"},
    {"ApplyV3ToVector", ApplyV3ToVector, METH_VARARGS, "Apply the v3 to a given vector"},
    {"ApplyV4ToDyn", ApplyV4ToDyn, METH_VARARGS, "Apply the v3 to a given dynamical matrix"},
    {"ApplyV3_FT", ApplyV3_FT, METH_VARARGS, "Apply the full v3 at finite temperature"},
    {"ApplyV4_FT", ApplyV4_FT, METH_VARARGS, "Apply the full v4 at finite temperature"},
    {NULL, NULL, 0, NULL}
};


// Module initialization
PyMODINIT_FUNC initsscha_HP_odd(void) {
    (void) Py_InitModule("sscha_HP_odd", odd_engine);
}

// --------------------------------------------------------
/*
 * Here the functions of the module.
 */

/*
 * APPLY D3 TO DYN (VECTOR)
 * ===============
 * 
 *  Apply the D3 to a vector that represents a dynamical matrix
 * This function returns an output vector.
 * 
 * This function is accelerated using OpenMP
 * 
 * Parameters
 * ----------
 *  X : double vector (size = (n_modes, n_random))
 *  Y : double vector (size = (n_modes, n_random))
 *  rho : double vector (size = nrandom)
 *  w : double vector (size = n_modes)
 *  T : double 
 *  input_dyn : double vector (size = n_modes^2)  
 *  output_vector : double vector (size = n_modes)
 *  mode : int 
 *    If mode = 1 : OpenMP version
 *    If mode = 2 : MPI version (Not yet implemented)
 */
static PyObject *ApplyV3ToDyn(PyObject * self, PyObject * args) {
  PyArrayObject * npy_X, *npy_Y, *npy_rho, *npy_omega, *npy_input, *npy_output, *npy_symmetries, *npy_n_deg, *npy_deg_space;
  double * X; 
  double *Y;
  double *w;
  double *input;
  double *output;
  double * symmetries;
  int * n_deg;
  int ** good_deg_space;
  double * rho;
  int N_configs;
  int N_modes;
  int mode;
  double T;

  int index_mode = 0, index_config = 1;

  // Parse the python arguments
  if (!PyArg_ParseTuple(args, "OOOOdOOiOOO", &npy_X, &npy_Y, &npy_rho, &npy_omega, &T, &npy_input, &npy_output, &mode,
			&npy_symmetries, &npy_n_deg, &npy_deg_space))
    return NULL;
  
  // Check the array memory setting
  if (npy_X->flags & NPY_ARRAY_F_CONTIGUOUS) {
    index_mode = 1;
    index_config = 0;
  }
  

  // Get the dimension of the arrays
  N_modes = PyArray_DIM(npy_X, index_mode);
  N_configs = PyArray_DIM(npy_X, index_config);

  // Check the dimensions of all the variables
  if (N_configs != PyArray_DIM(npy_rho,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_configs from X is %d, while len(rho) = %d\n", N_configs, PyArray_DIM(npy_rho, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_omega,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_modes from X is %d, while len(w) = %d\n", N_modes, PyArray_DIM(npy_omega, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_output,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes, PyArray_DIM(npy_output, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes*N_modes != PyArray_DIM(npy_input,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes*N_modes, PyArray_DIM(npy_input, 0));
    exit(EXIT_FAILURE);
  }

  // Retrive the pointer to the data from the python object
  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rho = (double*) PyArray_DATA(npy_rho);
  w = (double*) PyArray_DATA(npy_omega);
  input = (double*) PyArray_DATA(npy_input);
  output = (double*) PyArray_DATA(npy_output);

  // Read the symmetries
  symmetries = (double*)PyArray_DATA(npy_symmetries);
  n_deg = (int*)PyArray_DATA(npy_n_deg);

  // Build the degeneracy space
  good_deg_space = (int **) malloc(sizeof(int*) * N_modes);
  int i, j;
  int counter= 0;
  int N_symmetries;
  for (i = 0; i < N_modes;++i) {
    good_deg_space[i] = (int*) malloc(sizeof(int) * n_deg[i]);
    for (j = 0; j < n_deg[i]; ++j) {
      good_deg_space[i][j] = ((int*) PyArray_DATA(npy_deg_space))[counter++];
    }
  }

  N_symmetries = PyArray_DIM(npy_symmetries, 0);


  // Check the mode
  if (mode == 1) {
    OMP_ApplyD3ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  } else if (mode == 2) {
    // Use the MPI version
    MPI_ApplyD3ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  }
  else {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "mode %d not implemented.\n", mode);
    exit(EXIT_FAILURE);
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *ApplyV3_FT(PyObject * self, PyObject * args) {
  PyArrayObject * npy_X, *npy_Y, *npy_rho, *npy_omega, *npy_input, *npy_output, *npy_symmetries, *npy_n_deg, *npy_deg_space;
  double * X; 
  double *Y;
  double *w;
  double *input;
  double *output;
  double * symmetries;
  int * n_deg;
  int ** good_deg_space;
  double * rho;
  int N_configs;
  int N_modes;
  int mode, start_A, end_A;
  int transpose;
  double T;

  int index_mode = 0, index_config = 1;

  // Parse the python arguments
  if (!PyArg_ParseTuple(args, "OOOOdOOiOOOiii", &npy_X, &npy_Y, &npy_rho, &npy_omega, &T, &npy_input, &npy_output, &mode,
			&npy_symmetries, &npy_n_deg, &npy_deg_space, &start_A, &end_A, &transpose))
    return NULL;
  
  // Check the array memory setting
  if (npy_X->flags & NPY_ARRAY_F_CONTIGUOUS) {
    index_mode = 1;
    index_config = 0;
  }
  

  // Get the dimension of the arrays
  N_modes = PyArray_DIM(npy_X, index_mode);
  N_configs = PyArray_DIM(npy_X, index_config);

  // Check the dimensions of all the variables
  if (N_configs != PyArray_DIM(npy_rho,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_configs from X is %d, while len(rho) = %d\n", N_configs, PyArray_DIM(npy_rho, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_omega,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_modes from X is %d, while len(w) = %d\n", N_modes, PyArray_DIM(npy_omega, 0));
    exit(EXIT_FAILURE);
  }
  if (end_A != PyArray_DIM(npy_output,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", end_A, PyArray_DIM(npy_output, 0));
    exit(EXIT_FAILURE);
  }
  if (end_A != PyArray_DIM(npy_input,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", end_A, PyArray_DIM(npy_input, 0));
    exit(EXIT_FAILURE);
  }

  // Retrive the pointer to the data from the python object
  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rho = (double*) PyArray_DATA(npy_rho);
  w = (double*) PyArray_DATA(npy_omega);
  input = (double*) PyArray_DATA(npy_input);
  output = (double*) PyArray_DATA(npy_output);

  // Read the symmetries
  symmetries = (double*)PyArray_DATA(npy_symmetries);
  n_deg = (int*)PyArray_DATA(npy_n_deg);

  // Build the degeneracy space
  good_deg_space = (int **) malloc(sizeof(int*) * N_modes);
  int i, j;
  int counter= 0;
  int N_symmetries;

  printf("Degenerate space in C:\n");
  for (i = 0; i < N_modes;++i) {
    good_deg_space[i] = (int*) malloc(sizeof(int) * n_deg[i]);
      printf("Mode %d -> ", i);
    for (j = 0; j < n_deg[i]; ++j) {
      good_deg_space[i][j] = ((int*) PyArray_DATA(npy_deg_space))[counter++];
      printf(" %d ", good_deg_space[i][j]);
    }
    printf("\n");
  }
  fflush(stdout);

  N_symmetries = PyArray_DIM(npy_symmetries, 0);


  // Check the mode
  //if (mode == 1) {
    //OMP_ApplyD3ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  //  exit(EXIT_FAILURE);
  if (mode == 2) {
    // Use the MPI version
    MPI_D3_FT(X, Y, rho, w, T, N_modes, start_A, end_A, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space,
      transpose);
  }
  else {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "mode %d not implemented.\n", mode);
    exit(EXIT_FAILURE);
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *ApplyV4_FT(PyObject * self, PyObject * args) {
  PyArrayObject * npy_X, *npy_Y, *npy_rho, *npy_omega, *npy_input, *npy_output, *npy_symmetries, *npy_n_deg, *npy_deg_space;
  double * X; 
  double *Y;
  double *w;
  double *input;
  double *output;
  double * symmetries;
  int * n_deg;
  int ** good_deg_space;
  double * rho;
  int N_configs;
  int N_modes;
  int mode, start_A, end_A;
  double T;
  int transpose;

  int index_mode = 0, index_config = 1;

  // Parse the python arguments
  if (!PyArg_ParseTuple(args, "OOOOdOOiOOOiii", &npy_X, &npy_Y, &npy_rho, &npy_omega, &T, &npy_input, &npy_output, &mode,
			&npy_symmetries, &npy_n_deg, &npy_deg_space, &start_A, &end_A, &transpose))
    return NULL;
  
  // Check the array memory setting
  if (npy_X->flags & NPY_ARRAY_F_CONTIGUOUS) {
    index_mode = 1;
    index_config = 0;
  }
  

  // Get the dimension of the arrays
  N_modes = PyArray_DIM(npy_X, index_mode);
  N_configs = PyArray_DIM(npy_X, index_config);

  // Check the dimensions of all the variables
  if (N_configs != PyArray_DIM(npy_rho,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_configs from X is %d, while len(rho) = %d\n", N_configs, PyArray_DIM(npy_rho, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_omega,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_modes from X is %d, while len(w) = %d\n", N_modes, PyArray_DIM(npy_omega, 0));
    exit(EXIT_FAILURE);
  }
  if (end_A != PyArray_DIM(npy_output,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", end_A, PyArray_DIM(npy_output, 0));
    exit(EXIT_FAILURE);
  }
  if (end_A != PyArray_DIM(npy_input,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", end_A, PyArray_DIM(npy_input, 0));
    exit(EXIT_FAILURE);
  }

  // Retrive the pointer to the data from the python object
  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rho = (double*) PyArray_DATA(npy_rho);
  w = (double*) PyArray_DATA(npy_omega);
  input = (double*) PyArray_DATA(npy_input);
  output = (double*) PyArray_DATA(npy_output);

  // Read the symmetries
  symmetries = (double*)PyArray_DATA(npy_symmetries);
  n_deg = (int*)PyArray_DATA(npy_n_deg);

  // Build the degeneracy space
  good_deg_space = (int **) malloc(sizeof(int*) * N_modes);
  int i, j;
  int counter= 0;
  int N_symmetries;
  for (i = 0; i < N_modes;++i) {
    good_deg_space[i] = (int*) malloc(sizeof(int) * n_deg[i]);
    for (j = 0; j < n_deg[i]; ++j) {
      good_deg_space[i][j] = ((int*) PyArray_DATA(npy_deg_space))[counter++];
    }
  }

  N_symmetries = PyArray_DIM(npy_symmetries, 0);


  // Check the mode
  //if (mode == 1) {
    //OMP_ApplyD3ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  //  exit(EXIT_FAILURE);
  if (mode == 2) {
    // Use the MPI version
    MPI_D4_FT(X, Y, rho, w, T, N_modes, start_A, end_A, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space, transpose);
  }
  else {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "mode %d not implemented.\n", mode);
    exit(EXIT_FAILURE);
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject *ApplyV3ToVector(PyObject * self, PyObject * args) {
  PyArrayObject * npy_X, *npy_Y, *npy_rho, *npy_omega, *npy_input, *npy_output, *npy_symmetries, *npy_n_deg, *npy_deg_space;
  double * X, *Y, *w, *input, *output, *rho;
  double * symmetries;
  int * n_deg;
  double * deg_space;
  int ** good_deg_space;
  int N_configs, N_modes, mode;
  double T;

  int index_mode = 0, index_config = 1;

  // Parse the python arguments
  if (!PyArg_ParseTuple(args, "OOOOdOOiOOO", &npy_X, &npy_Y, &npy_rho, &npy_omega, &T, &npy_input, &npy_output, &mode,
			&npy_symmetries, &npy_n_deg, &npy_deg_space))
    return NULL;

  // Check the array memory setting
  if (npy_X->flags & NPY_ARRAY_F_CONTIGUOUS) {
    index_mode = 1;
    index_config = 0;
  }
  
  // Get the dimension of the arrays
  N_modes = PyArray_DIM(npy_X, index_mode);
  N_configs = PyArray_DIM(npy_X, index_config);

  // Check the dimensions of all the variables
  if (N_configs != PyArray_DIM(npy_rho,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_configs from X is %d, while len(rho) = %d\n", N_configs, PyArray_DIM(npy_rho, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_omega,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_modes from X is %d, while len(w) = %d\n", N_modes, PyArray_DIM(npy_omega, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes*N_modes != PyArray_DIM(npy_output,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes*N_modes, PyArray_DIM(npy_output, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_input,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes, PyArray_DIM(npy_input, 0));
    exit(EXIT_FAILURE);
  }

  // Retrive the pointer to the data from the python object
  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rho = (double*) PyArray_DATA(npy_rho);
  w = (double*) PyArray_DATA(npy_omega);
  input = (double*) PyArray_DATA(npy_input);
  output = (double*) PyArray_DATA(npy_output);

  
  // Read the symmetries
  symmetries = (double*)PyArray_DATA(npy_symmetries);
  n_deg = (int*)PyArray_DATA(npy_n_deg);

  // Build the degeneracy space
  good_deg_space = (int**) malloc(sizeof(int*) * N_modes);
  int i, j;
  int counter= 0;
  int N_symmetries;
  for (i = 0; i < N_modes;++i) {
    good_deg_space[i] = (int*) malloc(sizeof(int) * n_deg[i]);
    for (j = 0; j < n_deg[i]; ++j) {
      good_deg_space[i][j] = ((int*) PyArray_DATA(npy_deg_space))[counter++];
    }
  }

  N_symmetries = PyArray_DIM(npy_symmetries, 0);
  
  // Check the mode
  if (mode == 1) {
    //printf("I'm here\n");
    fflush(stdout);
    OMP_ApplyD3ToVector(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  } else if (mode == 2) {
    MPI_ApplyD3ToVector(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  }
  else {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "mode %d not implemented.\n", mode);
    exit(EXIT_FAILURE);
  }

  Py_INCREF(Py_None);
  return Py_None;
}



/*
 * APPLY D4 TO DYN 
 * ===============
 * 
 *  Apply the D4 to a vector that represents a dynamical matrix
 * This function returns an output vector.
 * 
 * This function is accelerated using OpenMP
 * 
 * Parameters
 * ----------
 *  X : double vector (size = (n_modes, n_random))
 *  Y : double vector (size = (n_modes, n_random))
 *  rho : double vector (size = nrandom)
 *  w : double vector (size = n_modes)
 *  T : double 
 *  input_dyn : double vector (size = n_modes^2)  
 *  output_vector : double vector (size = n_modes)
 *  mode : int 
 *    If mode = 1 : OpenMP version
 *    If mode = 2 : MPI version (Not yet implemented)
 */
static PyObject *ApplyV4ToDyn(PyObject * self, PyObject * args) {
  PyArrayObject * npy_X, *npy_Y, *npy_rho, *npy_omega, *npy_input, *npy_output, *npy_symmetries, *npy_n_deg, *npy_deg_space;
  double * X; 
  double *Y;
  double *w;
  double *input;
  double *output;
  double * rho;
  int N_configs;
  int N_modes;
  int mode;
  double T;
  double * symmetries;
  int * n_deg;
  double * deg_space;
  int ** good_deg_space;


  int index_mode = 0, index_config = 1;

  // Parse the python arguments
  if (!PyArg_ParseTuple(args, "OOOOdOOiOOO", &npy_X, &npy_Y, &npy_rho, &npy_omega, &T, &npy_input, &npy_output, &mode,
			&npy_symmetries, &npy_n_deg, &npy_deg_space))
    return NULL;
  
  // Check the array memory setting
  if (npy_X->flags & NPY_ARRAY_F_CONTIGUOUS) {
    index_mode = 1;
    index_config = 0;
  }
  

  // Get the dimension of the arrays
  N_modes = PyArray_DIM(npy_X, index_mode);
  N_configs = PyArray_DIM(npy_X, index_config);

  // Check the dimensions of all the variables
  if (N_configs != PyArray_DIM(npy_rho,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_configs from X is %d, while len(rho) = %d\n", N_configs, PyArray_DIM(npy_rho, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes != PyArray_DIM(npy_omega,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "N_modes from X is %d, while len(w) = %d\n", N_modes, PyArray_DIM(npy_omega, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes*N_modes != PyArray_DIM(npy_output,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes*N_modes, PyArray_DIM(npy_output, 0));
    exit(EXIT_FAILURE);
  }
  if (N_modes*N_modes != PyArray_DIM(npy_input,0)) {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "The output vector should have a length of %d instead of %d\n", N_modes*N_modes, PyArray_DIM(npy_input, 0));
    exit(EXIT_FAILURE);
  }

  // Retrive the pointer to the data from the python object
  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rho = (double*) PyArray_DATA(npy_rho);
  w = (double*) PyArray_DATA(npy_omega);
  input = (double*) PyArray_DATA(npy_input);
  output = (double*) PyArray_DATA(npy_output);

    
  // Read the symmetries
  symmetries = (double*)PyArray_DATA(npy_symmetries);
  n_deg = (int*)PyArray_DATA(npy_n_deg);

  // Build the degeneracy space
  good_deg_space = (int**) malloc(sizeof(int*) * N_modes);
  int i, j;
  int counter= 0;
  int N_symmetries;
  for (i = 0; i < N_modes;++i) {
    good_deg_space[i] = (int*) malloc(sizeof(int) * n_deg[i]);
    for (j = 0; j < n_deg[i]; ++j) {
      good_deg_space[i][j] = ((int*) PyArray_DATA(npy_deg_space))[counter++];
    }
  }

  N_symmetries = PyArray_DIM(npy_symmetries, 0);


  // Check the mode
  if (mode == 1) {
    OMP_ApplyD4ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space);
  } else if (mode == 2) {
    MPI_ApplyD4ToDyn(X, Y, rho, w, T, N_modes, N_configs, input, output, symmetries, N_symmetries, n_deg, good_deg_space); 
  }
  else {
    fprintf(stderr, "Error in file %s, line %d:\n", __FILE__ ,  __LINE__);
    fprintf(stderr, "mode %d not implemented.\n", mode);
    exit(EXIT_FAILURE);
  }

  Py_INCREF(Py_None);
  return Py_None;
}



/*
 * GET V3
 * ======
 * 
 * This function computes the v3 matrix in real space by using the replica symmetry algorithm.
 * The v3 is compute in the basis of the modes.  
 * 
 *  X_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *  Y_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *
 * Parameters
 * ----------
 *      X : ndarray(size = (N_modes, N_random))
 *          The array of the position (in the mode basis)
 *      Y : ndarray(size = (N_modes, N_random))
 *          The array of the forces (in the mode basis)
 *      N_modes : int
 *          The first dimension of the X/Y arrays
 *      N_random : int
 *          The second dimension of the X/Y arrrays
 *      out_v3 : ndarray(size = (N_modes, N_modes, N_modes))
 *          The output array that will be filled with the v3
 */

// #ifdef _MPI
// static PyObject * GetV3(PyObject * self, PyObject * args) {
//   int N_modes, N_random;
//   PyObject * npy_X, *npy_Y, *npy_v3;
//   double *X;
//   double *Y;
//   double *v3;
//   double tmp1, tmp2;
  
//   // Get the path dir
//   if (!PyArg_ParseTuple(args, "OOiiO", &npy_X, &npy_Y, &N_modes, &N_random, &npy_v3))
//     return NULL;

//   // Convert the array
//   //numpy_array = PyArray_FROM_OTF(bare_array, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

//   X = (double*) PyArray_DATA(npy_X);
//   Y = (double*) PyArray_DATA(npy_Y);
//   v3 = (double*) PyArray_DATA(npy_v3);


//   int i;
//   int rank=0, size=1;

//   // Check the parallelization and assign the size and the rank
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);

//   int a, b, c;
//   for (a = 0; a < N_modes; ++a) {
//     for (b = 0; b < N_modes; ++b) {
//       for (c = 0; c < N_modes; ++c) {
//         // Compute the averages
//         tmp1 = 0;
//         for (i = rank; i < N_random; i += size) {
//           tmp1 += 
//             X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] +
//             X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] +
//             Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i];
//         }
//         tmp2 = tmp1;
//         // Perform the reduction
//         MPI_Allreduce(&tmp1, &tmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
//         v3[a * N_modes*N_modes + b * N_modes + c] = - tmp2 / (3*N_random); 
         
//       }
//     }
//   }


//   Py_INCREF(Py_None);
//   return Py_None;
// }
// #endif


/*
 * GET V3 (OPEN MP VERSION)
 * ======
 * 
 * This function computes the v3 matrix in real space by using the replica symmetry algorithm.
 * The v3 is compute in the basis of the modes.  
 * 
 *  X_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *  Y_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *
 * Parameters
 * ----------
 *      X : ndarray(size = (N_modes, N_random))
 *          The array of the position (in the mode basis)
 *      Y : ndarray(size = (N_modes, N_random))
 *          The array of the forces (in the mode basis)
 *      N_modes : int
 *          The first dimension of the X/Y arrays
 *      N_random : int
 *          The second dimension of the X/Y arrrays
 *      out_v3 : ndarray(size = (N_modes, N_modes, N_modes))
 *          The output array that will be filled with the v3
 */

//#ifndef _MPI
static PyObject * GetV3(PyObject * self, PyObject * args) {
  int N_modes, N_random;
  PyObject * npy_X, *npy_Y, *npy_v3;
  double *X;
  double *Y;
  double *v3;
  double tmp2;
  
  // Get the path dir
  if (!PyArg_ParseTuple(args, "OOiiO", &npy_X, &npy_Y, &N_modes, &N_random, &npy_v3))
    return NULL;

  // Convert the array
  //numpy_array = PyArray_FROM_OTF(bare_array, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  v3 = (double*) PyArray_DATA(npy_v3);


  int i;
  int rank=0, size=1;
  
  // Check the parallelization and assign the size and the rank

  int a, b, c;
#pragma omp parallel for collapse(3) shared(X, Y, v3) private(a,b,c)
  for (a = 0; a < N_modes; ++a) {
    for (b = 0; b < N_modes; ++b) {
      for (c = 0; c < N_modes; ++c) {
        // Compute the averages
	double tmp1 = 0;
	for (i = 1; i < N_random; i++) {
	  tmp1 += 
	    X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] +
	    X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] +
	    Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i];
	}

        v3[a * N_modes*N_modes + b * N_modes + c] = - tmp1 / (3*N_random); 
         
      }
    }
  }


  Py_INCREF(Py_None);
  return Py_None;
}
//#endif


/*
 * APPLY V4 (OpenMP version)
 * =========================
 * 
 * This function applies the supermatrix (1 - v4 Lambda) to
 * a vector
 * 
 *  X_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *  Y_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *  Lambda_\mu\nu
 *
 * Parameters
 * ----------
 *      X : ndarray(size = (N_modes, N_random), dtype = double)
 *          The array of the position (in the mode basis)
 *      Y : ndarray(size = (N_modes, N_random), dtype = double)
 *          The array of the forces (in the mode basis)
 *      Lambda : ndarray(size = (N_modes, N_modes), dtype = double)
 *          The Lambda matrix
 *      v_in : ndarray(size = (N_modes, N_modes), dtype = double)
 *          The vector to which the matrix is multiplied
 *      N_modes : int
 *          The first dimension of the X/Y arrays
 *      N_random : int
 *          The second dimension of the X/Y arrrays
 *      N_eff : double
 *          The sum over all the weights of the configurations
 *      out_v3 : ndarray(size = (N_modes, N_modes, N_modes))
 *          The output array that will be filled with the v3
 */

//#ifndef _MPI
static PyObject * ApplyV4(PyObject * self, PyObject * args) {
  int N_modes, N_random;
  PyObject * npy_X, *npy_Y, *npy_vin, *npy_vout;
  PyObject * npy_Lambda;
  double *X;
  double *Y;
  double * Lambda;
  double *v_in, *v_out;
  double N_eff;
  double v4;
  
  // Get the path dir
  if (!PyArg_ParseTuple(args, "OOOOiidO", &npy_X, &npy_Y, &npy_Lambda, &npy_vin, 
      &N_modes, &N_random, &N_eff, &npy_vout))
    return NULL;

  // Convert the array
  //numpy_array = PyArray_FROM_OTF(bare_array, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  Lambda = (double*) PyArray_DATA(npy_Lambda);
  v_in = (double*) PyArray_DATA(npy_vin);
  v_out = (double*) PyArray_DATA(npy_vout);

  int i;
  int a, b, c, d;
  for (a = 0; a < N_modes; ++a) {
    for (b = 0; b < N_modes; ++b) {
      v_out[N_modes * a + b] = v_in[N_modes * a + b];
      for (c = 0; c < N_modes; ++c) {
        for (d = 0; c < N_modes; ++d) {
          // Compute the v4
          v4 = 0;
          #pragma omp parallel 
          {
            double tmp1 = 0;
            #pragma omp for
            for (i = 1; i < N_random; i++) {
              tmp1 += 
                X[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i] * Y[N_random*d + i] +
                X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] * X[N_random*d + i] +
                X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] * X[N_random*d + i] +
                Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i] * X[N_random*d + i];
            }

            // Perform the reduction
            #pragma omp critical 
              v4 += tmp1;
          }
          v4 /= (4 * N_eff);
          v_out[N_modes * a + b] -= v4 * Lambda[N_modes * c + d] * v_in[N_modes * c + d];
        }
      }
    }
  }


  Py_INCREF(Py_None);
  return Py_None;
}
//#endif 



/*
 * GET V4 dot LAMBDA (OpenMP version)
 * ==================================
 * 
 * This function computes the matrix v4 matrix
 * a vector
 * 
 *  X_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 *  Y_\mu^i = \sum_x e_\mu^x u_x^i \sqrt{m_i} \omega_\mu
 * 
 * 
 *
 * Parameters
 * ----------
 *      X : ndarray(size = (N_modes, N_random), dtype = double)
 *          The array of the position (in the mode basis)
 *      Y : ndarray(size = (N_modes, N_random), dtype = double)
 *          The array of the forces (in the mode basis)
 *      N_modes : int
 *          The first dimension of the X/Y arrays
 *      N_random : int
 *          The second dimension of the X/Y arrrays
 *      N_eff : double
 *          The sum over all the weights of the configurations
 *      thr : double
 *          The threshold below which it is considered zero.
 *      out_rows : ndarray(size = N_modes*N_modes)
 *          The rows indices
 *      out_cols : ndarray(size = N_modes*N_modes)
 *          The cols indices
 *      out_v4 : ndarray(size = (N_modes, N_modes, N_modes, N_modes))
 *          The data for each element in rows cols.
 * 
 * Returns
 * -------
 *      n_tot : int
 *        The number of elements different from zero.
 */

//#ifndef _MPI
static PyObject * GetV4(PyObject * self, PyObject * args) {
  int N_modes, N_random;
  PyObject * npy_X, *npy_Y, *npy_rows, *npy_cols, *npy_v4out, *npy_Lambda;
  double *X;
  double *Y;
  double *v4_out;
  int * rows, *cols;
  double N_eff, thr;
  double v4;
  
  // Get the path dir
  if (!PyArg_ParseTuple(args, "OOiiddOOO", &npy_X, &npy_Y, 
      &N_modes, &N_random, &N_eff, &thr, &npy_rows, &npy_cols, &npy_v4out))
    return NULL;

  // Convert the array
  //numpy_array = PyArray_FROM_OTF(bare_array, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  X = (double*) PyArray_DATA(npy_X);
  Y = (double*) PyArray_DATA(npy_Y);
  rows = (int*) PyArray_DATA(npy_rows);
  cols = (int*) PyArray_DATA(npy_cols);
  v4_out = (double*) PyArray_DATA(npy_v4out);

  int i;
  int a, b, c, d;
  int n_tot = 0;
  int n_equal;
  #pragma omp parallel shared(X, Y, rows, cols, v4_out, n_tot) private(a,b,c,d,v4, i)
  {
    #pragma omp for
    for (a = 0; a < N_modes; ++a) {
      for (b = a; b < N_modes; ++b) {
        for (c = b; c < N_modes; ++c) {
          for (d = c; d < N_modes; ++d) {
            // Compute the v4
            v4 = 0;
            for (i = 1; i < N_random; i++) {
              v4 += 
                X[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i] * Y[N_random*d + i] +
                X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] * X[N_random*d + i] +
                X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] * X[N_random*d + i] +
                Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i] * X[N_random*d + i];
            }
            v4 /= (4 * N_eff);


            // Add the element to the sparse matrix (and all the symmetries)
            if (fabs(v4) > thr) {

              // Count how many times the element will be repeated
              if (a == b && b == c && c == d) {
                v4 /= 24;
              }
              else if((a == b && b == c) || (a == b && b == d) || 
                (a == c && c == d) || (b == c && c == d) ) {
                v4 /= 6;
              }
              else if(a == b || a == c || a == d || b==c || b==d || c==d) {
                v4 /= 2;
              }

              // Here the critical part
              // Where we update the solution.
              // Only 1 process at time can use it.
              #pragma omp critical 
              {
                rows[n_tot] = N_modes * a + b;
                cols[n_tot] = N_modes * c + d;
                v4_out[n_tot] = v4;
                n_tot++;

                rows[n_tot] = N_modes * b + a;
                cols[n_tot] = N_modes * c + d;
                v4_out[n_tot] = v4;
                n_tot++;
              
                rows[n_tot] = N_modes * c + a;
                cols[n_tot] = N_modes * b + d;
                v4_out[n_tot] = v4;
                n_tot++;
              
                rows[n_tot] = N_modes * a + c;
                cols[n_tot] = N_modes * b + d;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * b + c;
                cols[n_tot] = N_modes * a + d;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * c + b;
                cols[n_tot] = N_modes * a + d;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * c + b;
                cols[n_tot] = N_modes * d + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * b + c;
                cols[n_tot] = N_modes * d + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + c;
                cols[n_tot] = N_modes * b + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * c + d;
                cols[n_tot] = N_modes * b + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * b + d;
                cols[n_tot] = N_modes * c + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + b;
                cols[n_tot] = N_modes * c + a;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + a;
                cols[n_tot] = N_modes * c + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * a + d;
                cols[n_tot] = N_modes * c + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * c + d;
                cols[n_tot] = N_modes * a + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + c;
                cols[n_tot] = N_modes * a + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * a + c;
                cols[n_tot] = N_modes * d + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * c + a;
                cols[n_tot] = N_modes * d + b;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * b + a;
                cols[n_tot] = N_modes * d + c;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * a + b;
                cols[n_tot] = N_modes * d + c;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + b;
                cols[n_tot] = N_modes * a + c;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * b + d;
                cols[n_tot] = N_modes * a + c;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * a + d;
                cols[n_tot] = N_modes * b + c;
                v4_out[n_tot] = v4;
                n_tot++;
                
                rows[n_tot] = N_modes * d + a;
                cols[n_tot] = N_modes * b + c;
                v4_out[n_tot] = v4;
                n_tot++;
              }
            }
          }
        }
      }
    }
  }

  return PyInt_FromLong(n_tot);
}
//#endif 
