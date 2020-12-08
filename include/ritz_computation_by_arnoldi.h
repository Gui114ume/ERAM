#ifndef RITZ_BY_ARNOLDI_H
#define RITZ_BY_ARNOLDI_H 

#include <lapacke.h>
#include "types.h"
#include "matrix_computation.h"

void eigenvalues_print(int order, double * wr, double * wi);

void eigenvector_print(char * vector_type, int order, double * wi, double * vector, int ldv);

void ritz_computation_by_arnoldi(Matrix M, double * wr, double * wi, double * vr);

double residu_ritz(Matrix M, double * wr, double * vr, int order);

#endif