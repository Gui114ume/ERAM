#ifndef ALGO_GRAMSHMIT_H
#define ALGO_GRAMSHMIT_H 

#include "types.h"
#include "matrix_computation.h"
#include "rdtsc.h"
#include "parameters.h"

Arnoldi_res Arnoldi_modified(Matrix A, Vector v, int n, int m);

#endif