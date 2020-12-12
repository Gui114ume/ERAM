#ifndef MATRIX_COMPUTATION_H
#define MATRIX_COMPUTATION_H 

#include <math.h>
#include "types.h"

void Matrix_print (Matrix M);

void Vector_print (Vector V);

Matrix dot_product_matrix(Matrix A, Matrix B);

void dot_prod_matrix_vector_from_matrix(Vector vout, Matrix A, Matrix B, int column);

double dot_product_vector(Vector A, Vector B);

double dot_product_vector_matrix_column(Vector A, Matrix B, int column);

Vector prod_matrix_column_coef(Matrix A, int column, double k);

void vector_minus_vector(Vector A, Vector B);

double norm_vector(Vector A);

double norm_matrix_column(Matrix A, int column);

void normalize_vector(Vector A);

int check_orthogolality(Matrix A);

void from2d_to_1d(double *M_oneD, Matrix M);

Matrix coeff_dot_identy(double k, int order);

Matrix matrix_minus_matrix(Matrix A, Matrix B);

Vector dot_Matrix_Vector(Matrix A, double* vr, int order, int index);

double vector_norm(Vector v);

void sort_eigen(double* wr, double* vr, int order);

void recompute_initial_vector_explicit(Vector initial, int order, double *vr);

void recompute_initial_vector_implicit(Matrix H, Matrix V, Ritz_eigen ritz, int shift_p);

void QRfactorisation(Matrix M, Matrix Q, Matrix R);

void vector_dot_coeff(Vector V, double k);

Matrix transposed(Matrix M);

void dot_prod_mat_vec(double * res , Matrix A, double * vect);

double scalaire_vec_vec(double * vout_1, double * vout_2, int n);

double norme_vecteur_double(double * vout, int n);
#endif
