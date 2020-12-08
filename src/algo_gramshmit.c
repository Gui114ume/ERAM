#include "algo_gramshmit.h"


Arnoldi_res Arnoldi_modified(Matrix A, Vector v, int n, int m)
{
	// DEFINE -------------------------------------------------
	if (n > m )
	{
		printf("Error : n > m , n = %d m = %d\n", n, m);								
		exit(1);
	}
	if (A.X_SIZE != A.Y_SIZE)
	{
		printf("Error : Entry matrix not : m*m\n");
		exit(1);

	}
	Arnoldi_res res;
	
	res.V.X_SIZE = n;
	res.V.Y_SIZE = m;
	res.V.A = malloc(sizeof(double*)*res.V.Y_SIZE);
	for (int i = 0; i < res.V.Y_SIZE; ++i)
	{
		res.V.A[i] = calloc(res.V.X_SIZE, sizeof(double));
	}
	
	res.H.X_SIZE = m;
	res.H.Y_SIZE = m;
	res.H.A = malloc(sizeof(double*)*res.H.Y_SIZE);
	for (int i = 0; i < res.H.Y_SIZE; ++i)
	{
		res.H.A[i] = calloc(res.H.X_SIZE, sizeof(double));
	}
	Vector vtmp;
	vtmp.Y_SIZE = A.Y_SIZE;
	vtmp.V = malloc(sizeof(double)*vtmp.Y_SIZE);
	// INITIALIZATION -----------------------------------------
	
	normalize_vector(v);
	
	for (int i = 0; i < res.V.Y_SIZE; ++i)
	{
		res.V.A[i][0] = v.V[i];
	}
	// ALGORITHM ----------------------------------------------
	for (int j = 0; j < n-1; ++j)
	{
		vtmp = dot_prod_matrix_vector_from_matrix(A, res.V, j);
		for (int i = 0; i < j+1; ++i)
		{
			res.H.A[i][j] = dot_product_vector_matrix_column(vtmp,res.V,i);
			vtmp = vector_minus_vector(vtmp,
									   prod_matrix_column_coef(res.V,i,res.H.A[i][j]));
		}
		res.H.A[j+1][j] = norm_vector(vtmp);
		
		for (int i = 0; i < res.V.Y_SIZE; ++i)
		{
			res.V.A[i][j+1] = vtmp.V[i]/res.H.A[j+1][j];
		}
	}
	
	return res;

}
