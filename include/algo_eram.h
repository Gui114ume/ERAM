#ifndef ALGO_ERAM_H
#define ALGO_ERAM_H

#include "types.h"
#include "matrix_computation.h"
#include "rdtsc.h"
#include "parameters.h"


int find_max(double* module,
	     int nb_values);
  
Eigenvalues n_biggest(Eigenvalues eigenvalues,
		      int nb);

Eigenvalues compute_eigenvalues(Matrix_1d hessenberg);

Eigenvectors compute_eigenvectors(Matrix_1d   hessenberg,
				  Eigenvalues eigenvalues);

Vector create_initial_vector(Eigenvectors eigenvectors);

int condition(double rs);

void Eram();

#endif
