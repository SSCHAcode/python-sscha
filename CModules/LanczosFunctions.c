#include "LanczosFunctions.h"

#define DEB 0

// The eigenvalues of the Covariance matrix
double f_ups(double w, double T) {
    double n_w = 0;
    if (T > 0) {
        n_w = 1 / ( exp(w * RY_TO_K / T) - 1);
    }
    return 2 *w / (1 + n_w);
}



void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes,
			 int N_configs, const double * input_vector, double * output_dyn,
			 double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space ) {

    // Compute the N_eff
    double N_eff = 0;
    int i;

    //#pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    if (DEB) {
      printf("File %s, Line %d: N_eff = %.4f\n", __FILE__, __LINE__, N_eff);
      fflush(stdout);
    }
    
    // Prepare the new modified X
    double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    for (i = 0; i < N_configs*N_modes; ++i) {
        new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    }

    if (DEB) {
      printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__, N_eff);
      fflush(stdout);
    }
    

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    // Perform the application
    int a, b, c, new_a, new_b, new_c;
    int j, k, i_sym, N_sym_tmp;
    double sym_coeff = 0;
    //#pragma omp parallel for collapse(3) private(a,b,c)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
	        // Check if this element is zero by symmetry
	        int stop= 0;
	        for (i = 0; i < N_sym; ++i) {
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		if (stop == 1) continue;
	      
                double tmp = 0;
		
		//#pragma omp parallel for private(i) reduce(+:tmp)
                for (i = 0; i < N_configs; ++i) {
                    tmp += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
                    //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
                    //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];

                }


		// Apply all the symmetries in the degenerate subspace
		for (i = 0; i < N_degeneracy[a]; ++i) {
		  new_a = degenerate_space[a][i];
		  for (j = 0; j < N_degeneracy[b]; ++j) {
		    new_b = degenerate_space[b][j];
		    for (k = 0; k < N_degeneracy[c]; ++k) {
		      new_c = degenerate_space[c][k];

		      // Check if there are degeneracies
		      // If not, symmetries are useless, apply only the identity
		      N_sym_tmp = N_sym;
		      if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] == 1)
			N_sym_tmp = 1;

		      for (i_sym = 0; i_sym < N_sym; ++i_sym) {
			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

			
			output_dyn[new_a * N_modes + new_b] += -tmp * input_vector[new_c] * sym_coeff / (3 * N_eff * N_sym);
			output_dyn[new_b * N_modes + new_a] += -tmp * input_vector[new_c] * sym_coeff / (3 * N_eff * N_sym);
			output_dyn[new_a * N_modes + new_c] += -tmp * input_vector[new_b] * sym_coeff / (3 * N_eff * N_sym);
			output_dyn[new_c * N_modes + new_a] += -tmp * input_vector[new_b] * sym_coeff / (3 * N_eff * N_sym);
			output_dyn[new_c * N_modes + new_b] += -tmp * input_vector[new_a] * sym_coeff / (3 * N_eff * N_sym);
			output_dyn[new_b * N_modes + new_c] += -tmp * input_vector[new_a] * sym_coeff / (3 * N_eff * N_sym);
		      }
		    }
		  }
		}

		// Check if the symmetries are mixing something
		
		if (DEB && c == N_modes -1) {
		  printf("a = %d, b = %d: output = %.8e\n", a, b, output_dyn[a * N_modes + b]);
		}
	    }
	   
        }
    }

    // Free memory
    free(new_X);
}



void OMP_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

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

		if (DEB && b == N_modes - 1 && c == N_modes - 1) {
		  printf("a = %d,  output = %.8e\n", a, output_vector[a]);
		}
	
            }
        }
    }

    // Free memory
    free(new_X);
}



void OMP_ApplyD4ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
    double T, int N_modes, int N_configs, const double * input_dyn, double * output_dyn){

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
    int a, b, c, d;
#pragma omp parallel for collapse(4) private(a,b,c, d)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
	      for (d = 0; d < N_modes; ++d) {
		
                double tmp = 0;
                for (i = 0; i < N_configs; ++i) {
		  double tmp1;
		  tmp1 = new_X[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i] * Y[N_configs*d +i];
		  tmp1 += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * new_X[N_configs*d +i];
		  tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  tmp += tmp1 * rho[i];
                }

                output_dyn[a*N_modes + b] += -tmp * input_dyn[N_modes * b + c] / (4 * N_eff);
	      }	
            }
        }
    }

    // Free memory
    free(new_X);
}



