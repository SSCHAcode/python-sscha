#ifndef LANCZOS_FUNCTION_H
#define LANCZOS_FUNCTION_H

#include <stdio.h>
#include <math.h>

#define RY_TO_K 157887.32400374097

/*
 * Apply the D3 matrix to a vector (OpenMP)
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
 */
void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, 
    double T, int N_modes, int N_configs, const double * input_vector, double * output_dyn);


/*
 * Apply the D3 matrix to a dyn (OpenMP)
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
 */
void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, 
    double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector);



#endif