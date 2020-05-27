#include "LanczosFunctions.h"


#define DEB 1

// The eigenvalues of the Covariance matrix
double f_ups(double w, double T) {
    double n_w = 0;
    if (T > 0) {
        n_w = 1 / ( exp(w * RY_TO_K / T) - 1);
    }
    return 2 *w / (1 + n_w);
}



// void OMP_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes,
// 			 int N_configs, const double * input_vector, double * output_dyn,
// 			 double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space ) {

//     // Compute the N_eff
//     double N_eff = 0;
//     int i;

//     //#pragma omp parallel for private(i) reduction(+:N_eff)
//     for (i = 0; i < N_configs; ++i)
//         N_eff += rho[i];

//     if (DEB) {
//       printf("File %s, Line %d: N_eff = %.4f\n", __FILE__, __LINE__, N_eff);
//       fflush(stdout);
//     }
    
//     // Prepare the new modified X
//     double * new_X = malloc(sizeof(double) * N_configs* N_modes);

//     //#pragma omp parallel for private(i)
//     for (i = 0; i < N_configs*N_modes; ++i) {
//         new_X[i] = X[i] * f_ups(w[i / N_configs], T);
//     }

//     if (DEB) {
//       printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__, N_eff);
//       fflush(stdout);
//     }
    

//     // Initialize the output
//     for (i = 0; i < N_modes*N_modes; ++i)
//         output_dyn[i] = 0;

//     if (DEB) printf("Applying the d3 to vector!\n");

//     // Perform the application
//     int a, b, c, new_a, new_b, new_c;
//     int j, k, i_sym, N_sym_tmp;
//     double sym_coeff = 0;
//     //#pragma omp parallel for collapse(3) private(a,b,c)
//     for (a = 0; a < N_modes; ++a) {
//         for (b = 0; b < N_modes; ++b) {
//             for (c = 0; c < N_modes; ++c) {
// 	        // Check if this element is zero by symmetry
// 	        int stop= 0;

// 		if (DEB) printf("Element a=%d, b=%d, c=%d ... \n", a, b, c);
// 	        for (i = 0; i < N_sym; ++i) {
// 		  if (DEB) printf("%d) s_aa = %.2f, s_bb = %.2f, s_cc = %.2f\n", i,
// 			 symmetries[i * N_modes*N_modes + a*N_modes + a],
// 			 symmetries[i * N_modes*N_modes + b*N_modes + b],
// 			 symmetries[i * N_modes*N_modes + c*N_modes + c]);
// 		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
// 			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
// 			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
// 		    stop = 1;
// 		    break;
// 		  }
// 		}
// 		if (stop == 1) continue;
// 		if (DEB)printf("I'm computing this element.\n");
	      
//                 double tmp = 0;
		
// 		//#pragma omp parallel for private(i) reduce(+:tmp)
//                 for (i = 0; i < N_configs; ++i) {
//                     tmp += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
//                     //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
//                     //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];

//                 }


// 		// Apply all the symmetries in the degenerate subspace
// 		for (i = 0; i < N_degeneracy[a]; ++i) {
// 		  new_a = degenerate_space[a][i];
// 		  for (j = 0; j < N_degeneracy[b]; ++j) {
// 		    new_b = degenerate_space[b][j];
// 		    for (k = 0; k < N_degeneracy[c]; ++k) {
// 		      new_c = degenerate_space[c][k];

// 		      // Check if there are degeneracies
// 		      // If not, symmetries are useless, apply only the identity
// 		      N_sym_tmp = N_sym;
// 		      if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] == 1)
// 			N_sym_tmp = 1;

// 		      if (DEB)
// 			printf("Deg space = %d | new_a = %d, new_b = %d, new_c = %d\n", N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c],
// 			       new_a, new_b, new_c);
		      
// 		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
// 			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
// 			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
// 			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

// 			if (DEB)
// 			  printf("IN_VEC_OUT_DYN: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));


			
// 			output_dyn[new_a * N_modes + new_b] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_dyn[new_b * N_modes + new_a] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_dyn[new_a * N_modes + new_c] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_dyn[new_c * N_modes + new_a] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_dyn[new_c * N_modes + new_b] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_dyn[new_b * N_modes + new_c] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 		      }
// 		    }
// 		  }
// 		}
// 		if (DEB)
// 		  printf("\n");

// 		// Check if the symmetries are mixing something
		
// 		if (DEB && c == N_modes -1) {
// 		  printf("a = %d, b = %d: output = %.8e\n", a, b, output_dyn[a * N_modes + b]);
// 		}
// 	    }
	   
//         }
//     }

//     // Free memory
//     free(new_X);
// }


/*
 * This function uses the new defintions of the X and Y matrices
 */
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
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    if (DEB) {
      printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__, N_eff);
      fflush(stdout);
    }
    

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    if (DEB) printf("Applying the d3 to vector!\n");

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

		if (DEB) printf("Element a=%d, b=%d, c=%d ... \n", a, b, c);
	        for (i = 0; i < N_sym; ++i) {
		  if (DEB) printf("%d) s_aa = %.2f, s_bb = %.2f, s_cc = %.2f\n", i,
			 symmetries[i * N_modes*N_modes + a*N_modes + a],
			 symmetries[i * N_modes*N_modes + b*N_modes + b],
			 symmetries[i * N_modes*N_modes + c*N_modes + c]);
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		if (stop == 1) continue;
		if (DEB)printf("I'm computing this element.\n");
	      
                double tmp = 0;
		
		//#pragma omp parallel for private(i) reduce(+:tmp)
                for (i = 0; i < N_configs; ++i) {
                    tmp += X[N_configs*a + i] * X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
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

		      if (DEB)
			printf("Deg space = %d | new_a = %d, new_b = %d, new_c = %d\n", N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c],
			       new_a, new_b, new_c);
		      
		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

			if (DEB)
			  printf("IN_VEC_OUT_DYN: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));


			
			output_dyn[new_a * N_modes + new_b] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_dyn[new_b * N_modes + new_a] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_dyn[new_a * N_modes + new_c] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_dyn[new_c * N_modes + new_a] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_dyn[new_c * N_modes + new_b] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_dyn[new_b * N_modes + new_c] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
		      }
		    }
		  }
		}
		if (DEB)
		  printf("\n");

		// Check if the symmetries are mixing something
		
		if (DEB && c == N_modes -1) {
		  printf("a = %d, b = %d: output = %.8e\n", a, b, output_dyn[a * N_modes + b]);
		}
	    }
	   
        }
    }

    // Free memory
    //free(new_X);
}



void MPI_ApplyD3ToVector(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes,
			 int N_configs, double * input_vector, double * output_dyn,
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
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    if (DEB) {
      printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__);
      fflush(stdout);
    }

	// Allocate a new output (This is used to perform a reduction of the MPI processors)
	double * new_output = (double*) calloc(sizeof(double),  N_modes * N_modes);
    

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    if (DEB) printf("Applying the d3 to vector!\n");

    // Perform the application
    int a, b, c, new_a, new_b, new_c;
    int j, k, i_sym, N_sym_tmp;
    double sym_coeff = 0;

	// MPI parallelization
	// NOTE MPI must be initialized
	int size=1, rank=0;
	int count, remainer, start, stop;
	#ifdef _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the data to everyone
	MPI_Bcast(input_vector, N_modes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif


	// The workload for each MPI process
	count = (N_modes*N_modes*N_modes) / size;
	// If the number of MPI process does not match, distribute them correctly
	remainer = (N_modes*N_modes*N_modes) % size;

	// Distribute the work in a clever way 
	if (rank < remainer) {
		start = rank * (count + 1);
		stop = start + count + 1;
	} else {
		start = rank * count + remainer;
		stop = start + count;
	}

	if (DEB) printf("Division: %d in %d ranks, with %d remainer\n",
		count*size, size, remainer);

	int mpi_index;
	for (mpi_index = start; mpi_index < stop; ++mpi_index) {
		c = mpi_index % N_modes;
		b = (mpi_index/ N_modes) % N_modes;
		a = (mpi_index/N_modes) / N_modes;

		if (DEB) {
			printf("RANK %d, index = %d (%d, %d). a = %d, b = %d, c = %d\n",
				rank, mpi_index, start, stop, a, b, c);
		}

		// Check if this element is zero by symmetry
		int stop= 0;

		if (DEB) printf("Element a=%d, b=%d, c=%d ... \n", a, b, c);
	        for (i = 0; i < N_sym; ++i) {
		  if (DEB) printf("%d) s_aa = %.2f, s_bb = %.2f, s_cc = %.2f\n", i,
			 symmetries[i * N_modes*N_modes + a*N_modes + a],
			 symmetries[i * N_modes*N_modes + b*N_modes + b],
			 symmetries[i * N_modes*N_modes + c*N_modes + c]);
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		if (stop == 1) continue;
		if (DEB)printf("I'm computing this element.\n");
	
		double tmp = 0;
		
		for (i = 0; i < N_configs; ++i) {
			tmp += X[N_configs*a + i] * X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
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

		      if (DEB)
				printf("Deg space = %d | new_a = %d, new_b = %d, new_c = %d\n", N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c],
					new_a, new_b, new_c);
		      
		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
				sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  	symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  	symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

				if (DEB)
				printf("IN_VEC_OUT_DYN: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));


				
				new_output[new_a * N_modes + new_b] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
				new_output[new_b * N_modes + new_a] += -tmp * input_vector[new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
				new_output[new_a * N_modes + new_c] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
				new_output[new_c * N_modes + new_a] += -tmp * input_vector[new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
				new_output[new_c * N_modes + new_b] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
				new_output[new_b * N_modes + new_c] += -tmp * input_vector[new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
		      }
		    }
		  }
		}
	}

	// Reduce the output dyn
	#ifdef _MPI
	MPI_Allreduce(new_output, output_dyn, N_modes*N_modes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	#ifndef _MPI
	// Copy the new output inside the output file.
	for (i = 0; i < N_modes*N_modes; ++i)
		output_dyn[i] = new_output[i];
	#endif
	   
    // Free memory
    free(new_output);
}


// Apply the full D3 at finite temperature
void MPI_D3_FT(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes, int start_A, int end_A,
			 int N_configs, double * input_psi, double * output_psi,
			 double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space, int transpose) {

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
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    if (DEB) {
      printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__);
      fflush(stdout);
    }

	// Allocate a new output (This is used to perform a reduction of the MPI processors)
	double * new_output = (double*) calloc(sizeof(double),  end_A);
    

    // Initialize the output
    for (i = 0; i < end_A; ++i)
        output_psi[i] = 0;

    if (DEB) printf("Applying the d3 to vector!\n");

    // Perform the application
    int a, b, c, new_a, new_b, new_c;
    int j, k, i_sym, N_sym_tmp;
    double sym_coeff = 0;
	double * n_w = (double*) calloc(sizeof(double), N_modes);

	// Fill the boson occupation numbers
	if (T > 0) 
		for (i = 0; i < N_modes; ++i) 
			n_w[i] = 1.0 / (exp(w[i] / (K_B * T)) - 1);

	// MPI parallelization
	// NOTE MPI must be initialized
	int size=1, rank=0;
	int count, remainer;
	unsigned long long int start, stop;
	#ifdef _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the data to everyone
	MPI_Bcast(input_psi, end_A, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif





	// The workload for each MPI process
	count = (N_modes*N_modes*N_modes) / size;
	// If the number of MPI process does not match, distribute them correctly
	remainer = (N_modes*N_modes*N_modes) % size;

	// Distribute the work in a clever way 
	if (rank < remainer) {
		start = rank * (count + 1);
		stop = start + count + 1;
	} else {
		start = rank * count + remainer;
		stop = start + count;
	}

	if (DEB) printf("Division: %d in %d ranks, with %d remainer\n",
		count*size, size, remainer);

	printf("Fast D3 FT Computation | rank %d computes from %u to %u\n", rank, start, stop);
	fflush(stdout);

	clock_t d3_timing = 0, sym_timing = 0;
	clock_t tmp_timing;
	double mult_coeff;


	unsigned long long int mpi_index;
	for (mpi_index = start; mpi_index < stop; ++mpi_index) {
		c = mpi_index % N_modes;
		b = (mpi_index/ N_modes) % N_modes;
		a = (mpi_index/N_modes) / N_modes;

		if (DEB) {
			printf("RANK %d, index = %d (%d, %d). a = %d, b = %d, c = %d\n",
				rank, mpi_index, start, stop, a, b, c);
		}

		// Check if this element is zero by symmetry
		int stop= 0;


		if (DEB) printf("Element a=%d, b=%d, c=%d ... \n", a, b, c);
	        for (i = 0; i < N_sym; ++i) {
		  if (DEB) printf("%d) s_aa = %.2f, s_bb = %.2f, s_cc = %.2f\n", i,
			 symmetries[i * N_modes*N_modes + a*N_modes + a],
			 symmetries[i * N_modes*N_modes + b*N_modes + b],
			 symmetries[i * N_modes*N_modes + c*N_modes + c]);
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		if (stop == 1) continue;
		if (DEB)printf("I'm computing this element.\n");
	

		// This is the time consuming part!!!
		// N_degeneracy is usually below 10, while N_configs can be hundreds of thounsands
		double tmp = 0;
		
		tmp_timing = clock();
		for (i = 0; i < N_configs; ++i) {
			tmp += X[N_configs*a + i] * X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
		}
		d3_timing += clock() - tmp_timing;


		// Apply all the symmetries in the degenerate subspace
		tmp_timing = clock();
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

		      if (DEB)
				printf("Deg space = %d | new_a = %d, new_b = %d, new_c = %d\n", N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c],
					new_a, new_b, new_c);
		      
		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
				sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  	symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  	symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

				if (DEB)
				printf("IN_VEC_OUT_DYN: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));

				// Here we must apply all the terms that contain d3
				// Get the final d3 with the correct symmetry coefficient
				tmp = -tmp * sym_coeff /  (3 * N_eff * N_sym_tmp);

				// We start applying them on R to fill the Y values of the output
				// This must take into account the transposition
				if (transpose == 0) mult_coeff = Z_coeff(w[new_a], n_w[new_a], w[new_b], n_w[new_b]);
				else mult_coeff = X2_coeff(w[new_a], n_w[new_a], w[new_b], n_w[new_b]);
				new_output[index_Y(new_a, new_b, N_modes)] += tmp * input_psi[new_c] * mult_coeff;

				if (DEB)
				printf("L_OP[ %d; %d] = %e\n", index_Y(new_a, new_b, N_modes), new_c, tmp * mult_coeff);

				if (transpose == 0) mult_coeff = Z_coeff(w[new_a], n_w[new_a], w[new_c], n_w[new_c]);
				else mult_coeff = X2_coeff(w[new_a], n_w[new_a], w[new_c], n_w[new_c]);
				new_output[index_Y(new_a, new_c, N_modes)] += tmp * input_psi[new_b] * mult_coeff;

				if (DEB)
				printf("L_OP[ %d; %d] = %e\n", index_Y(new_a, new_c, N_modes), new_b, tmp * mult_coeff);

				if (transpose == 0) mult_coeff = Z_coeff(w[new_c], n_w[new_c], w[new_b], n_w[new_b]);
				else  mult_coeff = X2_coeff(w[new_c], n_w[new_c], w[new_b], n_w[new_b]);
				new_output[index_Y(new_b, new_c, N_modes)] += tmp * input_psi[new_a] * mult_coeff;

				if (DEB)
				printf("L_OP[ %d; %d] = %e\n", index_Y(new_b, new_c, N_modes), new_a, tmp * mult_coeff);

				// We now apply on R to fill the A values of the output (Y2 is zero)
				if (transpose == 0) {
					new_output[index_A(new_a, new_b, N_modes)] += tmp * input_psi[new_c] * Z1_coeff(w[new_a], n_w[new_a], w[new_b], n_w[new_b]);
					new_output[index_A(new_a, new_c, N_modes)] += tmp * input_psi[new_b] * Z1_coeff(w[new_a], n_w[new_a], w[new_c], n_w[new_c]);
					new_output[index_A(new_b, new_c, N_modes)] += tmp * input_psi[new_a] * Z1_coeff(w[new_c], n_w[new_c], w[new_b], n_w[new_b]);
				}

				// Finally we apply on Y to fill the R values of the output
				// We must pay attention to double counting frequencies with exchange of axis.
				int extra_count = 1;
				if (new_b != new_c) extra_count = 2;
				if (transpose == 0) mult_coeff = X2_coeff(w[new_b], n_w[new_b], w[new_c], n_w[new_c]);
				else mult_coeff = Z_coeff(w[new_b], n_w[new_b], w[new_c], n_w[new_c]);
				new_output[new_a] += tmp * input_psi[index_Y(new_b, new_c, N_modes)] * mult_coeff * extra_count;

				extra_count = 1;
				if (new_a != new_c) extra_count = 2;
				if (transpose == 0) mult_coeff =  X2_coeff(w[new_a], n_w[new_a], w[new_c], n_w[new_c]);
				else mult_coeff = Z_coeff(w[new_a], n_w[new_a], w[new_c], n_w[new_c]);
				new_output[new_b] += tmp * input_psi[index_Y(new_a, new_c, N_modes)] * mult_coeff * extra_count;

				extra_count = 1;
				if (new_a != new_b) extra_count = 2;
				if (transpose == 0) mult_coeff = X2_coeff(w[new_b], n_w[new_b], w[new_a], n_w[new_a]);
				else mult_coeff = Z_coeff(w[new_b], n_w[new_b], w[new_a], n_w[new_a]);
				new_output[new_c] += tmp * input_psi[index_Y(new_b, new_a, N_modes)] * mult_coeff * extra_count;

		      }
		    }
		  }
		}
		sym_timing += clock() -tmp_timing;
	}

	printf("Fast D3 FT | rank %d | D3 timing: %.4lf sec | Sym timing: %.4lf sec\n", rank, d3_timing * 1.0 / CLOCKS_PER_SEC, sym_timing * 1.0 / CLOCKS_PER_SEC);
	fflush(stdout);

	// Reduce the output dyn
	#ifdef _MPI
	MPI_Allreduce(new_output, output_psi, end_A, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	#ifndef _MPI
	// Copy the new output inside the output file.
	for (i = 0; i < end_A; ++i)
		output_psi[i] = new_output[i];
	#endif
	   
    // Free memory
    free(new_output);
}



// Apply the full D3 at finite temperature
void MPI_D4_FT(const double * X, const double * Y, const double * rho, const double * w, double T, int N_modes, int start_A, int end_A,
			 int N_configs, double * input_psi, double * output_psi,
			 double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space, int transpose ) {

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
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    if (DEB) {
      printf("File %s, Line %d: Got the new X\n", __FILE__, __LINE__);
      fflush(stdout);
    }

	// Allocate a new output (This is used to perform a reduction of the MPI processors)
	double * new_output = (double*) calloc(sizeof(double),  end_A);
    

    // Initialize the output
    for (i = 0; i < end_A; ++i)
        output_psi[i] = 0;

    if (DEB) printf("Applying the d3 to vector!\n");

    // Perform the application
    int a, b, c, d, new_a, new_b, new_c, new_d;
    int j, k, h, i_sym, N_sym_tmp;
    double sym_coeff = 0;
	double * n_w = (double*) calloc(sizeof(double), N_modes);

	// Fill the boson occupation numbers
	if (T > 0) 
		for (i = 0; i < N_modes; ++i) 
			n_w[i] = 1.0 / (exp(w[i] / (K_B * T)) - 1);

	// MPI parallelization
	// NOTE MPI must be initialized
	int size=1, rank=0;
	unsigned long long int count, remainer, start, stop;
	#ifdef _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the data to everyone
	MPI_Bcast(input_psi, end_A, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif


	clock_t d4_timing = 0, sym_timing = 0;
	clock_t tmp_timing;

	int index1, index2;

	// The workload for each MPI process
	count = (N_modes*N_modes* (unsigned long long int) N_modes * N_modes) / size;
	// If the number of MPI process does not match, distribute them correctly
	remainer = (N_modes*N_modes* (unsigned long long int) N_modes * N_modes) % size;

	// Distribute the work in a clever way 
	if (rank < remainer) {
		start = rank * (count + 1);
		stop = start + count + 1;
	} else {
		start = rank * count + remainer;
		stop = start + count;
	}


	printf("Fast D4 FT Computation | rank %d computes from %u to %u\n", rank, start, stop);
	fflush(stdout);

	unsigned long long int mpi_index;
	for (mpi_index = start; mpi_index < stop; ++mpi_index) {
		
		d = mpi_index % N_modes;
		c = (mpi_index/N_modes) % N_modes;
		b = (mpi_index/(N_modes * N_modes)) % N_modes;
		a = (mpi_index/(N_modes*N_modes)) / N_modes;


		// Check if this element is zero by symmetry
		int stop= 0;

		for (i = 0; i < N_sym; ++i) {
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] *
			   symmetries[i * N_modes*N_modes + d*N_modes + d] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}

		if (stop == 1) continue;
		if (DEB)printf("I'm computing this element.\n");
	

		// This is the time consuming part!!!
		// N_degeneracy is usually below 10, while N_configs can be hundreds of thounsands
		double tmp = 0;
		tmp_timing = clock();
		for (i = 0; i < N_configs; ++i) {
			tmp += X[N_configs*a + i] * X[N_configs*b + i] * X[N_configs*c +i] * Y[N_configs*d +i] * rho[i];
		}
		d4_timing += clock() - tmp_timing;
		tmp_timing = clock();

		// Apply all the symmetries in the degenerate subspace
		for (i = 0; i < N_degeneracy[a]; ++i) {
		  new_a = degenerate_space[a][i];
		  for (j = 0; j < N_degeneracy[b]; ++j) {
		    new_b = degenerate_space[b][j];
		    for (k = 0; k < N_degeneracy[c]; ++k) {
		      new_c = degenerate_space[c][k];
			  for (h = 0; h < N_degeneracy[d]; ++h) {
				new_d = degenerate_space[d][h];

				// Check if there are degeneracies
				// If not, symmetries are useless, apply only the identity
				N_sym_tmp = N_sym;
				if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] * N_degeneracy[d] == 1)
					N_sym_tmp = 1;

				
				for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
					sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
						symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
						symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c] *
						symmetries[i_sym * N_modes * N_modes + d * N_modes + new_d];

					// Here we must apply all the terms that contain d3
					// Get the final d3 with the correct symmetry coefficient
					tmp = -tmp * sym_coeff /  (24 * N_eff * N_sym_tmp);

					int extra_count = 1; // (ab) (cd)
					if (new_c != new_d) extra_count = 2;
					new_output[index_Y(new_a, new_b, N_modes)] += tmp * input_psi[index_Y(new_c, new_d, N_modes)] * X_coeff(w[new_a], n_w[new_a],
						w[new_b], n_w[new_b], w[new_c], n_w[new_c], w[new_d], n_w[new_d]) * extra_count;
					
					extra_count = 1; // (ac) (bd)
					if (new_b != new_d) extra_count = 2;
					new_output[index_Y(new_a, new_c, N_modes)] += tmp * input_psi[index_Y(new_b, new_d, N_modes)] * X_coeff(w[new_a], n_w[new_a],
						w[new_c], n_w[new_c], w[new_b], n_w[new_b], w[new_d], n_w[new_d]) * extra_count;
					
					extra_count = 1; // (ad) (bc)
					if (new_b != new_c) extra_count = 2;
					new_output[index_Y(new_a, new_d, N_modes)] += tmp * input_psi[index_Y(new_b, new_c, N_modes)] * 
						X_coeff(w[new_a], n_w[new_a], w[new_d], n_w[new_d], w[new_b], n_w[new_b], w[new_c], n_w[new_c]) * extra_count;
					
					extra_count = 1; // (bc) (ad)
					if (new_a != new_d) extra_count = 2;
					new_output[index_Y(new_b, new_c, N_modes)] += tmp * input_psi[index_Y(new_a, new_d, N_modes)] * X_coeff(w[new_b], n_w[new_b],
						w[new_c], n_w[new_c], w[new_a], n_w[new_a], w[new_d], n_w[new_d]) * extra_count;

					extra_count = 1; // (bd) (ac)
					if (new_a != new_c) extra_count = 2;
					new_output[index_Y(new_b, new_d, N_modes)] += tmp * input_psi[index_Y(new_a, new_c, N_modes)] * X_coeff(w[new_b], n_w[new_b],
						w[new_d], n_w[new_d], w[new_a], n_w[new_a], w[new_c], n_w[new_c]) * extra_count;

					extra_count = 1; // (cd) (ab)
					if (new_a != new_b) extra_count = 2;
					new_output[index_Y(new_c, new_d, N_modes)] += tmp * input_psi[index_Y(new_a, new_b, N_modes)] * X_coeff(w[new_c], n_w[new_c],
						w[new_d], n_w[new_d], w[new_a], n_w[new_a], w[new_b], n_w[new_b]) * extra_count;



					// Now apply the X1
					extra_count = 1; // (ab) (cd)
					if (new_c != new_d) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_a, new_b, N_modes); index2 = index_Y(new_c, new_d, N_modes);}
					else {index2 = index_A(new_a, new_b, N_modes); index1 = index_Y(new_c, new_d, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * X1_coeff(w[new_a], n_w[new_a],
						w[new_b], n_w[new_b], w[new_c], n_w[new_c], w[new_d], n_w[new_d]) * extra_count;
					
					extra_count = 1; // (ac) (bd)
					if (new_b != new_d) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_a, new_c, N_modes); index2 = index_Y(new_b, new_d, N_modes);}
					else {index2 = index_A(new_a, new_c, N_modes); index1 = index_Y(new_b, new_d, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * X1_coeff(w[new_a], n_w[new_a],
						w[new_c], n_w[new_c], w[new_b], n_w[new_b], w[new_d], n_w[new_d]) * extra_count;
					
					extra_count = 1; // (ad) (bc)
					if (new_b != new_c) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_a, new_d, N_modes); index2 = index_Y(new_b, new_c, N_modes);}
					else {index2 = index_A(new_a, new_d, N_modes); index1 = index_Y(new_b, new_c, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * 
						X_coeff(w[new_a], n_w[new_a], w[new_d], n_w[new_d], w[new_b], n_w[new_b], w[new_c], n_w[new_c]) * extra_count;
					
					extra_count = 1; // (bc) (ad)
					if (new_a != new_d) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_b, new_c, N_modes); index2 = index_Y(new_a, new_d, N_modes);}
					else {index2 = index_A(new_b, new_c, N_modes); index1 = index_Y(new_a, new_d, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * X1_coeff(w[new_b], n_w[new_b],
						w[new_c], n_w[new_c], w[new_a], n_w[new_a], w[new_d], n_w[new_d]) * extra_count;

					extra_count = 1; // (bd) (ac)
					if (new_a != new_c) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_b, new_d, N_modes); index2 = index_Y(new_a, new_c, N_modes);}
					else {index2 = index_A(new_b, new_d, N_modes); index1 = index_Y(new_a, new_c, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * X1_coeff(w[new_b], n_w[new_b],
						w[new_d], n_w[new_d], w[new_a], n_w[new_a], w[new_c], n_w[new_c]) * extra_count;

					extra_count = 1; // (cd) (ab)
					if (new_a != new_b) extra_count = 2;
					if (transpose == 0) {index1 = index_A(new_c, new_d, N_modes); index2 = index_Y(new_a, new_b, N_modes);}
					else {index2 = index_A(new_c, new_d, N_modes); index1 = index_Y(new_a, new_b, N_modes);}
					new_output[index1] += tmp * input_psi[index2] * X1_coeff(w[new_c], n_w[new_c],
						w[new_d], n_w[new_d], w[new_a], n_w[new_a], w[new_b], n_w[new_b]) * extra_count;

				}
		      }
		    }
		  }
		}
		sym_timing += clock() - tmp_timing;
	}



	printf("Fast D4 FT | rank %d | D3 timing: %.4lf sec | Sym timing: %.4lf sec\n", rank, d4_timing * 1.0 / CLOCKS_PER_SEC, sym_timing * 1.0 / CLOCKS_PER_SEC);
	fflush(stdout);

	// Reduce the output dyn
	#ifdef _MPI
	MPI_Allreduce(new_output, output_psi, end_A, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	#ifndef _MPI
	// Copy the new output inside the output file.
	for (i = 0; i < end_A; ++i)
		output_psi[i] = new_output[i];
	#endif
	   
    // Free memory
    free(new_output);
}



// void OMP_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
// 		      double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector,
// 		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

//     // Compute the N_eff
//     double N_eff = 0;
//     int i;

//     //#pragma omp parallel for private(i) reduction(+:N_eff)
//     for (i = 0; i < N_configs; ++i)
//         N_eff += rho[i];

    
//     // Prepare the new modified X
//     double * new_X = malloc(sizeof(double) * N_configs* N_modes);

//     //#pragma omp parallel for private(i)
//     for (i = 0; i < N_configs*N_modes; ++i) {
//         new_X[i] = X[i] * f_ups(w[i / N_configs], T);
//     }

//     // Initialize the output
//     for (i = 0; i < N_modes; ++i)
//         output_vector[i] = 0;

//     // Perform the application
//     int a, b, c, new_a, new_b, new_c;
//     int j, k, i_sym, N_sym_tmp;
//     double sym_coeff = 0;
//     //#pragma omp parallel for collapse(3) private(a,b,c)
//     for (a = 0; a < N_modes; ++a) {
//         for (b = 0; b < N_modes; ++b) {
//             for (c = 0; c < N_modes; ++c) {
//                 double tmp = 0;

		
// 	        // Check if this element is zero by symmetry
// 	        int stop= 0;
// 	        for (i = 0; i < N_sym; ++i) {
// 		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
// 			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
// 			   symmetries[i * N_modes*N_modes + c*N_modes + c] + 1) < __EPSILON__) {
// 		    stop = 1;
// 		    break;
// 		  }
// 		}
		
// 		if (stop == 1) continue;
		
//                 for (i = 0; i < N_configs; ++i) {		 
//                     tmp += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
//                     //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
//                     //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];
//                     //tmp += tmp1 * rho[i];
//                 }
		
// 		// Apply all the symmetries in the degenerate subspace
// 		for (i = 0; i < N_degeneracy[a]; ++i) {
// 		  new_a = degenerate_space[a][i];
// 		  for (j = 0; j < N_degeneracy[b]; ++j) {
// 		    new_b = degenerate_space[b][j];
// 		    for (k = 0; k < N_degeneracy[c]; ++k) {
// 		      new_c = degenerate_space[c][k];

// 		      // Check if there are degeneracies
// 		      // If not, symmetries are useless, apply only the identity
// 		      N_sym_tmp = N_sym;
// 		      if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] == 1)
// 			N_sym_tmp = 1;

// 		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
// 			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
// 			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
// 			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

			
// 			if (DEB)
// 			  printf("IN_DYN_OUT_VEC: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));
			
			
// 			output_vector[new_a] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_vector[new_a] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_vector[new_b] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_vector[new_b] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_vector[new_c] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
// 			output_vector[new_c] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);

// 		      }
// 		    }
// 		  }
// 		}


// 		if (DEB && b == N_modes - 1 && c == N_modes - 1) {
// 		  printf("a = %d,  output = %.8e\n", a, output_vector[a]);
// 		}
	
//             }
//         }
//     }

//     // Free memory
//     free(new_X);
// }

void OMP_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs, const double * input_dyn, double * output_vector,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

    // Compute the N_eff
    double N_eff = 0;
    int i;

    //#pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    // Initialize the output
    for (i = 0; i < N_modes; ++i)
        output_vector[i] = 0;

    // Perform the application
    int a, b, c, new_a, new_b, new_c;
    int j, k, i_sym, N_sym_tmp;
    double sym_coeff = 0;
    //#pragma omp parallel for collapse(3) private(a,b,c)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
                double tmp = 0;

		
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
		
                for (i = 0; i < N_configs; ++i) {		 
                    tmp += X[N_configs*a + i] * X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
                    //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
                    //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];
                    //tmp += tmp1 * rho[i];
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

		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

			
			if (DEB)
			  printf("IN_DYN_OUT_VEC: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));
			
			
			output_vector[new_a] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_vector[new_a] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_vector[new_b] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_vector[new_b] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_vector[new_c] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			output_vector[new_c] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);

		      }
		    }
		  }
		}


		if (DEB && b == N_modes - 1 && c == N_modes - 1) {
		  printf("a = %d,  output = %.8e\n", a, output_vector[a]);
		}
	
            }
        }
    }

    // Free memory
    //free(new_X);
}


void MPI_ApplyD3ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs, double * input_dyn, double * output_vector,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

    // Compute the N_eff
    double N_eff = 0;
    int i;

    //#pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

	double * new_output = (double*)calloc(sizeof(double), N_modes);

    // Initialize the output
    for (i = 0; i < N_modes; ++i)
        output_vector[i] = 0;

    // Perform the application
    int a, b, c, new_a, new_b, new_c;
    int j, k, i_sym, N_sym_tmp;
    double sym_coeff = 0;


	// MPI parallelization
	// NOTE MPI must be initialized
	int size=1, rank=0;
	int count, remainer, start, stop;
	#ifdef _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the data to everyone
	MPI_Bcast(input_dyn, N_modes*N_modes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif


	// The workload for each MPI process
	count = (N_modes*N_modes*N_modes) / size;
	// If the number of MPI process does not match, distribute them correctly
	remainer = (N_modes*N_modes*N_modes) % size;

	// Distribute the work in a clever way 
	if (rank < remainer) {
		start = rank * (count + 1);
		stop = start + count + 1;
	} else {
		start = rank * count + remainer;
		stop = start + count;
	}

	int mpi_index;
	for (mpi_index = start; mpi_index < stop; ++mpi_index) {
		c = mpi_index % N_modes;
		b = (mpi_index/ N_modes) % N_modes;
		a = (mpi_index/N_modes) / N_modes;

		double tmp = 0;

		
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
		
                for (i = 0; i < N_configs; ++i) {		 
                    tmp += X[N_configs*a + i] * X[N_configs*b + i] * Y[N_configs*c +i] * rho[i];
                    //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i];
                    //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i];
                    //tmp += tmp1 * rho[i];
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

		      for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
			sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			  symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			  symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c];

			
			if (DEB)
			  printf("IN_DYN_OUT_VEC: symfactor = %.2f | d3[%d, %d, %d] = %.6e\n", sym_coeff, a, b, c, -tmp / (N_eff));
			
			
			new_output[new_a] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			new_output[new_a] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			new_output[new_b] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);
			new_output[new_b] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (6 * N_eff * N_sym_tmp);
			new_output[new_c] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (6 * N_eff * N_sym_tmp);
			new_output[new_c] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (6 * N_eff * N_sym_tmp);

		      }
		    }
		  }
		}


		if (DEB && b == N_modes - 1 && c == N_modes - 1) {
		  printf("a = %d,  output = %.8e\n", a, output_vector[a]);
		}
	
    }

	// Reduce the output dyn
	#ifdef _MPI
	MPI_Allreduce(new_output, output_vector, N_modes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	#ifndef _MPI
	// Copy the new output inside the output file.
	for (i = 0; i < N_modes; ++i)
		output_vector[i] = new_output[i];
	#endif
	   
    // Free memory
    free(new_output);
}




void OMP_ApplyD4ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs, const double * input_dyn, double * output_dyn,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

    // Compute the N_eff
    double N_eff = 0;
    int i;

    #pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    // Perform the application
    int a, b, c, d;
    int new_a, new_b, new_c, new_d;
    int j, k, h, i_sym, N_sym_tmp;
    double sym_coeff = 0;
    
    //#pragma omp parallel for collapse(4) private(a,b,c, d)
    for (a = 0; a < N_modes; ++a) {
        for (b = 0; b < N_modes; ++b) {
            for (c = 0; c < N_modes; ++c) {
	      for (d = 0; d < N_modes; ++d) {

		// Check if the operation is allowed by symmetries
		int stop = 0;
		
	        for (i = 0; i < N_sym; ++i) {
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] *
			   symmetries[i * N_modes*N_modes + d*N_modes + d] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		
		if (stop == 1) continue;
		
		
                double tmp = 0;
                for (i = 0; i < N_configs; ++i) {
		  tmp += X[N_configs*a + i] * X[N_configs*b + i] * X[N_configs*c +i] * Y[N_configs*d +i] * rho[i];
		  //tmp1 += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp += tmp1 * rho[i];
                }

		
		// Apply all the symmetries in the degenerate subspace
		for (i = 0; i < N_degeneracy[a]; ++i) {
		  new_a = degenerate_space[a][i];
		  for (j = 0; j < N_degeneracy[b]; ++j) {
		    new_b = degenerate_space[b][j];
		    for (k = 0; k < N_degeneracy[c]; ++k) {
		      new_c = degenerate_space[c][k];
		      for (h = 0; h < N_degeneracy[d]; ++h) {
			new_d = degenerate_space[d][h];

			// Check if there are degeneracies
			// If not, symmetries are useless, apply only the identity
			N_sym_tmp = N_sym;
			if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] * N_degeneracy[d] == 1)
			  N_sym_tmp = 1;

			for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
			  sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			    symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			    symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c] *
			    symmetries[i_sym * N_modes * N_modes + d * N_modes + new_d];

			
			if (DEB)
			  printf("IN_DYN_OUT_DYN: symfactor = %.2f | d4[%d, %d, %d, %d] = %.6e\n", sym_coeff, a, b, c, d, -tmp / (N_eff));
			
			
			output_dyn[N_modes * new_a + new_b] += -tmp * input_dyn[N_modes * new_c + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_b + new_a] += -tmp * input_dyn[N_modes * new_c + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_a + new_b] += -tmp * input_dyn[N_modes * new_d + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_b + new_a] += -tmp * input_dyn[N_modes * new_d + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);

			output_dyn[N_modes * new_a + new_c] += -tmp * input_dyn[N_modes * new_b + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_c + new_a] += -tmp * input_dyn[N_modes * new_b + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_a + new_c] += -tmp * input_dyn[N_modes * new_d + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_c + new_a] += -tmp * input_dyn[N_modes * new_d + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);

			output_dyn[N_modes * new_a + new_d] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_a] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_a + new_d] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_a] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);

			output_dyn[N_modes * new_b + new_c] += -tmp * input_dyn[N_modes * new_a + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_c + new_b] += -tmp * input_dyn[N_modes * new_a + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_b + new_c] += -tmp * input_dyn[N_modes * new_d + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_c + new_b] += -tmp * input_dyn[N_modes * new_d + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);

			output_dyn[N_modes * new_b + new_d] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_b] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_b + new_d] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_b] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);

			output_dyn[N_modes * new_c + new_d] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_c] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_c + new_d] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			output_dyn[N_modes * new_d + new_c] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			}
		      }
		    }
		  }
		}
		//output_dyn[a*N_modes + b] += -tmp * input_dyn[N_modes * b + c] / (4 * N_eff);
	      }	
            }
        }
    }

    // Free memory
    //free(new_X);
}



void MPI_ApplyD4ToDyn(const double * X, const double * Y, const double * rho, const double * w, 
		      double T, int N_modes, int N_configs,  double * input_dyn, double * output_dyn,
		      double * symmetries, int N_sym, int * N_degeneracy, int ** degenerate_space){

    // Compute the N_eff
    double N_eff = 0;
    int i;

    //#pragma omp parallel for private(i) reduction(+:N_eff)
    for (i = 0; i < N_configs; ++i)
        N_eff += rho[i];

    
    // Prepare the new modified X
    //double * new_X = malloc(sizeof(double) * N_configs* N_modes);

    //#pragma omp parallel for private(i)
    //for (i = 0; i < N_configs*N_modes; ++i) {
    //    new_X[i] = X[i] * f_ups(w[i / N_configs], T);
    //}

	double * new_output = (double*) calloc(sizeof(double), N_modes*N_modes);

    // Initialize the output
    for (i = 0; i < N_modes*N_modes; ++i)
        output_dyn[i] = 0;

    // Perform the application
    int a, b, c, d;
    int new_a, new_b, new_c, new_d;
    int j, k, h, i_sym, N_sym_tmp;
    double sym_coeff = 0;
    

	// MPI parallelization
	// NOTE MPI must be initialized
	int size=1, rank=0;
	long long count, start, stop;
	int remainer;
	#ifdef _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the data to everyone
	MPI_Bcast(input_dyn, N_modes*N_modes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//printf("MPI active\n");
    #endif
	#ifndef _MPI
	printf("Warning: MPI NOT ACTIVE! Check the compliation.\n");
	#endif


	// The workload for each MPI process
	count = (N_modes*N_modes * (long long) N_modes*N_modes) / size;
	// If the number of MPI process does not match, distribute them correctly
	remainer = (N_modes*N_modes* (long long) N_modes*N_modes) % size;

	// Distribute the work in a clever way 
	if (rank < remainer) {
		start = rank * (count + 1);
		stop = start + count + 1;
	} else {
		start = rank * count + remainer;
		stop = start + count;
	}

	// Print what we need to do for each processors
	printf("MPI process %d runs [%lld, %lld)\n", rank, start, stop);

	long long mpi_index;
	for (mpi_index = start; mpi_index < stop; ++mpi_index) {
		
		d = mpi_index % N_modes;
		c = (mpi_index/N_modes) % N_modes;
		b = (mpi_index/(N_modes * N_modes)) % N_modes;
		a = (mpi_index/(N_modes*N_modes)) / N_modes;

		// Check if the operation is allowed by symmetries
		int stop = 0;
		
	        for (i = 0; i < N_sym; ++i) {
		  if (fabs(symmetries[i * N_modes*N_modes + a*N_modes + a] *
			   symmetries[i * N_modes*N_modes + b*N_modes + b] *
			   symmetries[i * N_modes*N_modes + c*N_modes + c] *
			   symmetries[i * N_modes*N_modes + d*N_modes + d] + 1) < __EPSILON__) {
		    stop = 1;
		    break;
		  }
		}
		
		if (stop == 1) continue;
		
		
                double tmp = 0;
                for (i = 0; i < N_configs; ++i) {
		  tmp += X[N_configs*a + i] * X[N_configs*b + i] * X[N_configs*c +i] * Y[N_configs*d +i] * rho[i];
		  //tmp1 += new_X[N_configs*a + i] * new_X[N_configs*b + i] * Y[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp1 += new_X[N_configs*a + i] * Y[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp1 += Y[N_configs*a + i] * new_X[N_configs*b + i] * new_X[N_configs*c +i] * new_X[N_configs*d +i];
		  //tmp += tmp1 * rho[i];
                }

		
		// Apply all the symmetries in the degenerate subspace
		for (i = 0; i < N_degeneracy[a]; ++i) {
		  new_a = degenerate_space[a][i];
		  for (j = 0; j < N_degeneracy[b]; ++j) {
		    new_b = degenerate_space[b][j];
		    for (k = 0; k < N_degeneracy[c]; ++k) {
		      new_c = degenerate_space[c][k];
		      for (h = 0; h < N_degeneracy[d]; ++h) {
			new_d = degenerate_space[d][h];

			// Check if there are degeneracies
			// If not, symmetries are useless, apply only the identity
			N_sym_tmp = N_sym;
			if (N_degeneracy[a] * N_degeneracy[b] * N_degeneracy[c] * N_degeneracy[d] == 1)
			  N_sym_tmp = 1;

			for (i_sym = 0; i_sym < N_sym_tmp; ++i_sym) {
			  sym_coeff = symmetries[i_sym * N_modes * N_modes + a * N_modes + new_a] *
			    symmetries[i_sym * N_modes * N_modes + b * N_modes + new_b] *
			    symmetries[i_sym * N_modes * N_modes + c * N_modes + new_c] *
			    symmetries[i_sym * N_modes * N_modes + d * N_modes + new_d];

			
			if (DEB)
			  printf("IN_DYN_OUT_DYN: symfactor = %.2f | d4[%d, %d, %d, %d] = %.6e\n", sym_coeff, a, b, c, d, -tmp / (N_eff));
			
			
			new_output[N_modes * new_a + new_b] += -tmp * input_dyn[N_modes * new_c + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_b + new_a] += -tmp * input_dyn[N_modes * new_c + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_a + new_b] += -tmp * input_dyn[N_modes * new_d + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_b + new_a] += -tmp * input_dyn[N_modes * new_d + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);

			new_output[N_modes * new_a + new_c] += -tmp * input_dyn[N_modes * new_b + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_c + new_a] += -tmp * input_dyn[N_modes * new_b + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_a + new_c] += -tmp * input_dyn[N_modes * new_d + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_c + new_a] += -tmp * input_dyn[N_modes * new_d + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);

			new_output[N_modes * new_a + new_d] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_a] += -tmp * input_dyn[N_modes * new_b + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_a + new_d] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_a] += -tmp * input_dyn[N_modes * new_c + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);

			new_output[N_modes * new_b + new_c] += -tmp * input_dyn[N_modes * new_a + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_c + new_b] += -tmp * input_dyn[N_modes * new_a + new_d] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_b + new_c] += -tmp * input_dyn[N_modes * new_d + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_c + new_b] += -tmp * input_dyn[N_modes * new_d + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);

			new_output[N_modes * new_b + new_d] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_b] += -tmp * input_dyn[N_modes * new_a + new_c] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_b + new_d] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_b] += -tmp * input_dyn[N_modes * new_c + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);

			new_output[N_modes * new_c + new_d] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_c] += -tmp * input_dyn[N_modes * new_a + new_b] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_c + new_d] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			new_output[N_modes * new_d + new_c] += -tmp * input_dyn[N_modes * new_b + new_a] * sym_coeff / (24 * N_eff * N_sym_tmp);
			}
		      }
		    }
		  }
		}
		//output_dyn[a*N_modes + b] += -tmp * input_dyn[N_modes * b + c] / (4 * N_eff);
	      }	
	
	// Reduce the output dyn
	#ifdef _MPI
	MPI_Allreduce(new_output, output_dyn, N_modes*N_modes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif
	#ifndef _MPI
	// Copy the new output inside the output file.
	for (i = 0; i < N_modes*N_modes; ++i)
		output_dyn[i] = new_output[i];
	#endif
	   
    // Free memory
    free(new_output);

    // Free memory
    //free(new_X);
}


// The working functions
double Z_coeff(double w_a, double n_a, double w_b, double n_b) {
	return -2 * ((2*n_a + 1)*w_b + (2*n_b + 1)*w_a) / ((2*n_a + 1) * (2*n_b + 1));
}


double Z1_coeff(double w_a, double n_a, double w_b, double n_b) {
	return -2 * ( (2*n_a + 1)*w_b*n_b*(n_b + 1) + (2*n_b + 1)*w_a*n_a*(n_a+1)) / ((2*n_a + 1) * (2*n_b + 1));
}


double X2_coeff(double w_a, double n_a, double w_b, double n_b) {
	return -(2*n_b + 1) * (2*n_a +1) / (8*w_a *w_b);
}


double X_coeff(double w_a, double n_a, double w_b, double n_b, double w_c, double n_c, double w_d, double n_d) {
	return (2*n_c +1) * (2*n_d + 1)* (2*w_a * n_b + w_a + 2*w_b * n_a + w_b) / (4 * w_c * w_d * (2*n_a + 1)* (2*n_b + 1)); 
}

double X1_coeff(double w_a, double n_a, double w_b, double n_b, double w_c, double n_c, double w_d, double n_d) {
	return (2*n_c +1) * (2*n_d + 1)* (w_a*n_a*(n_a+1)*(2*n_b+1) + w_b*n_b*(n_b+1)*(2*n_a+1)) / (4 * w_c * w_d * (2*n_a + 1)* (2*n_b + 1)); 
}

int index_Y(int a, int b, int N) {
	// Get the symmetric index
	int ret = N;

	if (a >= b) ret += ((a - 1)*a) / 2 + b;
	else ret += ((b - 1)*b) / 2 + a;

	return ret;
}


int index_A(int a, int b, int N) {
	// Get the symmetric index
	int ret = N;

	// Go to Start-A
	ret += ((N+1)*N) / 2;

	if (a >= b) ret += ((a - 1)*a) / 2 + b;
	else ret += ((b - 1)*b) / 2 + a;

	return ret;
}