#include "LanczosFunctions.h"


// The eigenvalues of the Covariance matrix
double f_ups(double w, double T) {
    double n_w;
    if (T > 0) {
        n_w = 1 / ( exp(w * RY_TO_K / T) - 1);
    }
    return 2 *w / (1 + n_w);
}


void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes, int N_configs, const double * input_vector, double * output_dyn) {

    // Compute the N_eff
    double N_eff = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    #pragma omp parallel for private(i)
    for (i = 0; i < N_configs*N_modes; ++i) {
        new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    }

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    // Perform the application
    int a, b, c;
    #pragma omp parallel for collapse(3) private(a,b,c)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
                double tmp = 0;
                for (i = 0; i < N_configs; ++i) {
                    double tmp1;
                    tmp1 = new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i];
                    tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
                    tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];
                    tmp += tmp1 * rho[i];
                }

                output_dyn[a * N_modes + b] += -tmp * input_vector[c] / (3 * N_eff);
            }
        }
    }

    // Free memory
    free(new_X);
}



void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, 
    double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector){

    // Compute the N_eff
    double N_eff = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    #pragma omp parallel for private(i)
    for (i = 0; i < N_configs*N_modes; ++i) {
        new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    }

    // Initialize the output
    for (i = 0; i < N_modes; ++i)
        output_vector[i] = 0;

    // Perform the application
    int a, b, c;
    #pragma omp parallel for collapse(3) private(a,b,c)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
                double tmp = 0;
                for (i = 0; i < N_configs; ++i) {
                    double tmp1;
                    tmp1 = new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i];
                    tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
                    tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];
                    tmp += tmp1 * rho[i];
                }

                output_vector[a] += -tmp * input_dyn[N_modes * b + c] / (3 * N_eff);
            }
        }
    }

    // Free memory
    free(new_X);
}
    


