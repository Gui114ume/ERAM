#ifndef MATRIX_COMPUTATION_H
#define MATRIX_COMPUTATION_H 

#include "types.h"

void Matrix_print (Matrix M);

void Vector_print (Vector V);

Matrix dot_product_matrix(Matrix A, Matrix B);

Vector dot_prod_matrix_vector_from_matrix(Matrix A, Matrix B, int column);

double dot_product_vector(Vector A, Vector B);

double dot_product_vector_matrix_column(Vector A, Matrix B, int column);

Vector prod_matrix_column_coef(Matrix A, int column, double k);

Vector vector_minus_vector(Vector A, Vector B);

double norm_vector(Vector A);

void normalize_vector(Vector A);

int check_orthogolality(Matrix A);

#endif