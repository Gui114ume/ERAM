#ifndef ALGO_GRAMSHMIT_H
#define ALGO_GRAMSHMIT_H 

#include "types.h"
#include "matrix_computation.h"
#include "rdtsc.h"
#include "parameters.h"

void Arnoldi_modified(Arnoldi_res res, Matrix A, Vector v, int n, int m);

#endif