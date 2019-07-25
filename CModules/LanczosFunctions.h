#ifndef LANCZOS_FUNCTION_H
#define LANCZOS_FUNCTION_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RY_TO_K 157887.32400374097
#define __EPSILON__ 1e-6

//#define D_MPI
#ifdef _MPI
#include <mpi.h>
#endif

// NOTE: For now OpenMP parallelization is not active.

/*
 * Apply the D3 matrix to a vector (OpenMP/MPI)
 * ===============================
 * 
 * Apply the d3 matrix to a vector.
 * 
 * Both X and Y must have configurations as fast index
 * 
 * Parameters
 * ----------
 *      X : array (size (N_config * N_modes))
 *          The displacements of the ensemble in the mode representation
 *      Y : array (suze (N_config * N_modes))
 *          The forces of the ensemble in the mode representation
 *      rho : array( size = N_config)
 *          The weight of the configurations
 *      w : array(size = N_modes)
 *          The frequencies in Ry
 *      T : double
 *          The temperature of the calculation
 *      N_modes : int
 *      N_configs : int
 *      input_vector : array (size = N_modes)
 *          The input vector
 *      output_dyn : array (size = N_modes * N_modes)
 *          The result of the operation d3 * input_vector
 *      sym_matrix : array (size = (N_sym * N_modes * N_modes)
 *          The symmetry matrix in the mode basis.
 *      N_sym : int
 *          The number of symmetries
 *      N_degeneracy : array (size = N_modes)
 *          The degree of degeneracy of the mode
 *      degenerate_space : array of array [N_modes, N_degeneracy[i]]
 *          The mode indices that compose the degenerate subspace in which the first mode label belong to.
 */
void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, 
			 double T, int N_modes, int N_configs, const double * input_vector, double * output_dyn,
			 double * sym_matrix, int N_sym, int * N_degeneracy, int ** degenerate_space);
void MPI_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, 
			 double T, int N_modes, int N_configs,  double * input_vector, double * output_dyn,
			 double * sym_matrix, int N_sym, int * N_degeneracy, int ** degenerate_space);


/*
 * Apply the D3 matrix to a dyn (OpenMP/MPI)
 * ============================
 * 
 * Apply the d3 matrix to a dyn.
 * 
 * Parameters
 * ----------
 *      X : array (size (N_config * N_modes))
 *          The displacements of the ensemble in the mode representation
 *      Y : array (suze (N_config * N_modes))
 *          The forces of the ensemble in the mode representation
 *      rho : array( size = N_config)
 *          The weight of the configurations
 *      w : array(size = N_modes)
 *          The frequencies in Ry
 *      T : double
 *          The temperature of the calculation
 *      N_modes : int
 *      N_configs : int
 *      input_dyn : array (size = N_modes * N_modes)
 *          The input dyn
 *      output_vector : array (size = N_modes)
 *          The result of the operation d3 * input_dyn
 *      sym_matrix : array (size = (N_sym * N_modes * N_modes)
 *          The symmetry matrix in the mode basis.
 *      N_sym : int
 *          The number of symmetries
 *      N_degeneracy : array (size = N_modes)
 *          The degree of degeneracy of the mode
 *      degenerate_space : array of array [N_modes, N_degeneracy[i]]
 *          The mode indices that compose the degenerate subspace in which the first mode label belong to.
 */

void OMP_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
			 double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector,
			 double * sym_matrix, int N_sym, int * N_degeneracy, int ** degenerate_space);
void MPI_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
			 double T, int N_modes, int N_configs,  double * input_dyn, double * output_vector,
			 double * sym_matrix, int N_sym, int * N_degeneracy, int ** degenerate_space);



// Here the method to apply the d4 matrix (dyn to dyn)

void OMP_ApplyD4ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs, const double * input_dyn, double * output_dyn,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space);
void MPI_ApplyD4ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs,  double * input_dyn, double * output_dyn,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space);



#endif
