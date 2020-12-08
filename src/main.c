#include "types.h"
#include "algo_gramshmit.h"
#include "parameters.h"
#include "ritz_computation_by_arnoldi.h"

int main(int argc, char const *argv[])
{	srand(time(NULL));

	Arnoldi_res modified;
	Matrix M;
	Vector V;
	Ritz_eigen ritz;
	int order = 0;
	int iterations = 0;
	double * wr;
	double * wi;
	double * vr;

	double start_residu = 10000.;
	Vector residu; 
	residu.Y_SIZE = CONVERGENCE_ITERATIONS;
	residu.V = malloc(sizeof(double)*residu.Y_SIZE);
	//double convergence[CONVERGENCE_ITERATIONS];

	if (RANDOM_OR_VERIFIED_3X3)
	{
		M.X_SIZE = Matrix_size_X;
		M.Y_SIZE = Matrix_size_Y;
		V.Y_SIZE = Matrix_size_X;
		order = M.X_SIZE;
	}
	else
	{
		M.X_SIZE = 3;
		M.Y_SIZE = 3;
		V.Y_SIZE = 3;
		order = M.X_SIZE;
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

//Arnoldi modified --------------------------------------------

//Initialisation ----------------------------------------------
	printf("Modified Arnoldi\n");

	wr = malloc(sizeof(double)* order);
	wi = malloc(sizeof(double)* order);
	vr = malloc(sizeof(double)* order * order);

	ritz.order = order;

	ritz.eigen_value = malloc(sizeof(double) * ritz.order);
	ritz.eigen_vector = malloc(sizeof(Vector*) * ritz.order);
	for (int i = 0; i < ritz.order; ++i)
	{
		ritz.eigen_vector[i].Y_SIZE = ritz.order;
		ritz.eigen_vector[i].V = malloc(sizeof(double)* ritz.eigen_vector[i].Y_SIZE);
	}

//Algorithm run -----------------------------------------------
	if (EXPLICIT_RESTART)
	{
		while(start_residu > TOLERENCE && iterations < CONVERGENCE_ITERATIONS)
		{	
			modified = Arnoldi_modified(M,V,K_ITER,M.X_SIZE);
			
			ritz_computation_by_arnoldi(modified.H, wr, wi, vr);
			retype_sort_eigen(ritz, wr, wi, vr);
	
		   	start_residu = residu_ritz(M, ritz);
		   	residu.V[iterations] = start_residu;
		   	if (iterations > 1)
		   	{
		   		start_residu = abs(residu.V[iterations] - residu.V[iterations-1]);
		   		//convergence[iterations-1] = start_residu;
		   		//printf("Convergence %d : %f\n",iterations, convergence[iterations-1] );
		   	}
	
		   	recompute_initial_vector_explicit(V, ritz);
		   	iterations ++;
		}
	}

	//Print results -------------------------------------------
  	if (ALL_PRINT)
  	{	
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
		if (start_residu < TOLERENCE)
		{
			printf("Convergence YES, iterations %d \n", iterations);
	
		}
		else
			printf("Convergence NO, iterations %d\n", iterations);
  		//eigenvalues_print(order, wr, wi);
	   	//eigenvector_print("Right eigenvectors", order, wi, vr, order);
	   	//printf("Residu vector :\n");
	   	//Vector_print(residu);

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

	

	// FREES --------------------------------------------------
	for (int i = 0; i < M.Y_SIZE; i++)
	{
		free(M.A[i]);
	}
	free(M.A);
	free(V.V);

	return 0;
}