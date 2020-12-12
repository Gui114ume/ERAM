#include "algo_gramshmit.h"


void Arnoldi_modified(Arnoldi_res res, Matrix A, Vector v, int n, int m)
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

	
	res.V.X_SIZE = n;
	res.V.Y_SIZE = m;
	
	
	res.H.X_SIZE = m;
	res.H.Y_SIZE = m;
	
	Vector vtmp;
	Vector prod_matrix_vec;

	vtmp.Y_SIZE = A.Y_SIZE;
	vtmp.V = malloc(sizeof(double)*vtmp.Y_SIZE);
	// INITIALIZATION -----------------------------------------
	
	 double module = norm_vector(v);
	
	for (int i = 0; i < res.V.Y_SIZE; ++i)
	{
		res.V.A[i][0] = v.V[i]/module;
	}

	// ALGORITHM ----------------------------------------------
	for (int j = 0; j < n; ++j)
	{
		
		dot_prod_matrix_vector_from_matrix(vtmp, A, res.V, j);
	
		
		for (int i = 0; i < j+1; ++i)
		{
			res.H.A[i][j] = dot_product_vector_matrix_column(vtmp,res.V,i);
			
			prod_matrix_vec = prod_matrix_column_coef(res.V,i,res.H.A[i][j]);
			
			vector_minus_vector(vtmp, prod_matrix_vec);
							
			free(prod_matrix_vec.V);
		}
		res.H.A[j][j+1] = norm_vector(vtmp);
		if (res.H.A[j][j+1] == 0.0)
		{
			break;
		}
				
		if (j != n-1)
		{
			for (int i = 0; i < res.V.Y_SIZE; ++i)
			{
				res.V.A[i][j+1] = vtmp.V[i]/res.H.A[j][j+1];
			}
		}
	
		
	}
	free(vtmp.V);

	
}
