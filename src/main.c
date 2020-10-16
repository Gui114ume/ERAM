#include "types.h"
#include "algo_gramshmit.h"
#include "parameters.h"

int main(int argc, char const *argv[])
{	srand(time(NULL));

	Arnoldi_res classic;
	Arnoldi_res modified;
	Matrix M;
	Vector V;
	if (RANDOM_OR_VERIFIED_3X3)
	{
		M.X_SIZE = Matrix_size_X;
		M.Y_SIZE = Matrix_size_Y;
		V.Y_SIZE = Matrix_size_X;
	}
	else
	{
		M.X_SIZE = 3;
		M.Y_SIZE = 3;
		V.Y_SIZE = 3;
	}

	

	M.A = malloc(sizeof(double*)* M.Y_SIZE);
	V.V = malloc(sizeof(double)*V.Y_SIZE);

	for (int i = 0; i < M.Y_SIZE; i++)
	{
		M.A[i] = malloc(sizeof(double)*M.X_SIZE);
		
	}
if (RANDOM_OR_VERIFIED_3X3)
	{
				
		for (int i = 0; i < M.Y_SIZE; i++)
		{
			for (int j = 0; j < M.X_SIZE; ++j)
			{	
			
				M.A[i][j] = rand() % MAX_VALUE_MATRIX;
				
			}
		}
		for (int i = 0; i < V.Y_SIZE; ++i)
		{	
			V.V[i] = rand() % MAX_VALUE_MATRIX;		
		}

	}
	else
	{
		M.A[0][0] = 1;
		M.A[0][1] = 2;
		M.A[0][2] = 3;
		M.A[1][0] = 2;
		M.A[1][1] = 2;
		M.A[1][2] = 3;
		M.A[2][0] = 4;
		M.A[2][1] = 2;
		M.A[2][2] = 1;
		for (int i = 0; i < V.Y_SIZE; ++i)
		{	
			V.V[i] = 1.;	
		}
	}
if (INITIAL_MATRIX_AFF)
{
	printf("\nMatrix : \n");
	Matrix_print(M);
}
if (INITIAL_VECTOR_AFF)
{
	printf("Vector : \n");
	Vector_print(V);
}

//Arnoldi classic ---------------------------------------------
if (RUN_CLASSIC_MODIFIED == 0 || RUN_CLASSIC_MODIFIED == 2 )
{
	printf("Classic Arnoldi\n");
	classic = Arnoldi_classic(M,V,K_ITER,M.X_SIZE);
	if (H_MATRIX_AFF)
	{
		printf("Matrix H\n");
		Matrix_print(classic.H);
	}
	if (V_VECTOR_AFF)
	{
		printf("Vectors V\n");
		Matrix_print(classic.V);
	}
	if (ORTHOGONALITY_TEST)
	{	
		if (check_orthogolality(classic.V) == 0)
		{
			printf("Check orthogonality : YES\n\n" );
		}
		else
			printf("Check orthogonality : NO\n\n" );
	}
	// FREES --------------------------------------------------
	for (int i = 0; i < classic.H.Y_SIZE; ++i)
	{
		free(classic.H.A[i]);
	}
	free(classic.H.A);
	for (int i = 0; i < classic.V.Y_SIZE; ++i)
	{
		free(classic.V.A[i]);
	}
	free(classic.V.A);
}

//Arnoldi modified --------------------------------------------
if (RUN_CLASSIC_MODIFIED == 1 || RUN_CLASSIC_MODIFIED == 2 )
{
	printf("Modified Arnoldi\n");
	if (RANDOM_OR_VERIFIED_3X3)
	{
				
		for (int i = 0; i < M.Y_SIZE; i++)
		{
			for (int j = 0; j < M.X_SIZE; ++j)
			{	
			
				M.A[i][j] = rand() % MAX_VALUE_MATRIX;
				
			}
		}
		for (int i = 0; i < V.Y_SIZE; ++i)
		{	
			V.V[i] = rand() % MAX_VALUE_MATRIX;		
		}

	}
	else
	{
		M.A[0][0] = 1;
		M.A[0][1] = 2;
		M.A[0][2] = 3;
		M.A[1][0] = 2;
		M.A[1][1] = 2;
		M.A[1][2] = 3;
		M.A[2][0] = 4;
		M.A[2][1] = 2;
		M.A[2][2] = 1;
		for (int i = 0; i < V.Y_SIZE; ++i)
		{	
			V.V[i] = 1.;	
		}
	}
	
	modified = Arnoldi_modified(M,V,K_ITER,M.X_SIZE);
	if (H_MATRIX_AFF)
	{
		printf("Matrix H\n");
		Matrix_print(modified.H);
	}
	if (V_VECTOR_AFF)
	{
		printf("Vectors V\n");
		Matrix_print(modified.V);
	}
	if (ORTHOGONALITY_TEST)
	{	
		if (check_orthogolality(modified.V) == 0)
		{
			printf("Check orthogonality : YES\n\n" );
		}
		else
			printf("Check orthogonality : NO\n\n" );
	}

	// FREES --------------------------------------------------
	for (int i = 0; i < modified.H.Y_SIZE; ++i)
	{
		free(modified.H.A[i]);
	}
	free(modified.H.A);
	for (int i = 0; i < modified.V.Y_SIZE; ++i)
	{
		free(modified.V.A[i]);
	}
	free(modified.V.A);
}
	

	// FREES --------------------------------------------------
	for (int i = 0; i < M.Y_SIZE; i++)
	{
		free(M.A[i]);
	}
	free(M.A);
	free(V.V);

	return 0;
}