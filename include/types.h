#ifndef TYPES_H
#define TYPES_H 

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

struct matrix
{
	double **A;
	int X_SIZE;
	int Y_SIZE;
};
typedef struct matrix Matrix;

struct vector
{
	double *V;
	int Y_SIZE;
};
typedef struct vector Vector;

struct arnoldi_res
{
	Matrix V;
	Matrix H;
};
typedef struct arnoldi_res Arnoldi_res;

struct ritz_eigen
{
	int order;
	double * eigen_value;
	Vector * eigen_vector;
};
typedef struct ritz_eigen Ritz_eigen;

#endif