#include <python2.7/Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#ifdef _MPI
#include <mpi.h>
#endif

#include <omp.h>

static PyObject *GetV3 (PyObject * self, PyObject * args);


static PyMethodDef odd_engine[] = {
    {"GetV3", GetV3, METH_VARARGS, "Compute the v3 matrix using the replica method"},
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

#ifdef _MPI
static PyObject * GetV3(PyObject * self, PyObject * args) {
  int N_modes, N_random;
  PyObject * npy_X, *npy_Y, *npy_v3;
  double *X;
  double *Y;
  double *v3;
  double tmp1, tmp2;
  
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
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int a, b, c;
  for (a = 0; a < N_modes; ++a) {
    for (b = 0; b < N_modes; ++b) {
      for (c = 0; c < N_modes; ++c) {
        // Compute the averages
        tmp1 = 0;
        for (i = rank; i < N_random; i += size) {
          tmp1 += 
            X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] +
            X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] +
            Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i];
        }
        tmp2 = tmp1;
        // Perform the reduction
        MPI_Allreduce(&tmp1, &tmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        v3[a * N_modes*N_modes + b * N_modes + c] = - tmp2 / (3*N_random); 
         
      }
    }
  }


  Py_INCREF(Py_None);
  return Py_None;
}
#endif


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

#ifndef _MPI
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
  for (a = 0; a < N_modes; ++a) {
    for (b = 0; b < N_modes; ++b) {
      for (c = 0; c < N_modes; ++c) {
        // Compute the averages
        tmp2 = 0;
        #pragma omp parallel 
        {
          double tmp1 = 0;
          #pragma omp for
          for (i = 1; i < N_random; i++) {
            tmp1 += 
              X[N_random*a + i] * X[N_random*b + i] * Y[N_random*c + i] +
              X[N_random*a + i] * Y[N_random*b + i] * X[N_random*c + i] +
              Y[N_random*a + i] * X[N_random*b + i] * X[N_random*c + i];
          }

          // Perform the reduction
          #pragma omp critical 
            tmp2 += tmp1;
        }
        v3[a * N_modes*N_modes + b * N_modes + c] = - tmp2 / (3*N_random); 
         
      }
    }
  }


  Py_INCREF(Py_None);
  return Py_None;
}
#endif


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

#ifndef _MPI
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
        for (d = 0; c < N_modes; ++c) {
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
#endif 