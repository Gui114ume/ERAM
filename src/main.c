#include "types.h"
#include "algo_gramshmit.h"
#include "parameters.h"
#include "ritz_computation_by_arnoldi.h"

int main(int argc, char const *argv[])
{	
	srand(time(NULL));

	Matrix M;
	Vector V;
	int order = 0;
	int iterations = 0;
	double * wr;
	double * wi;
	double * vr;

	double start_residu = 10000.;
	Vector residu; 
	residu.Y_SIZE = CONVERGENCE_ITERATIONS;
	residu.V = malloc(sizeof(double)*residu.Y_SIZE);
	double convergence[CONVERGENCE_ITERATIONS];

	struct timeval timer_start, timer_stop; 

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
		M.A[0][1] = 3;
		M.A[0][2] = 3;
		M.A[1][0] = -2;
		M.A[1][1] = 11;
		M.A[1][2] = -2;
		M.A[2][0] = 8;
		M.A[2][1] = -7;
		M.A[2][2] = 6;
		for (int i = 0; i < V.Y_SIZE; ++i)
		{	
			V.V[i] = 1.;	
		}
	}
	if (ALL_PRINT)
	{
		printf("---------- ERAM ----------\n");
		printf("\n---------- ENTRY ---------\n");
		if (INITIAL_MATRIX_AFF)
		{
			printf("\nMatrix : n = %d ,m = %d\n",M.X_SIZE, M.Y_SIZE);
			Matrix_print(M);
		}
		if (INITIAL_VECTOR_AFF)
		{
			printf("Vector : \n");
			Vector_print(V);
		}
	}

//Arnoldi modified --------------------------------------------

//Initialisation ----------------------------------------------
	

	Arnoldi_res modified;
	
	wr = malloc(sizeof(double)* order);
	wi = malloc(sizeof(double)* order);
	vr = malloc(sizeof(double)* order * order);
/*
	Ritz_eigen ritz;

	ritz.order = order;

	ritz.eigen_value = (double*)malloc(sizeof(double) * order);
	ritz.eigen_vector = malloc(sizeof(Vector) * order);
	for (int i = 0; i < order; ++i)
	{
		ritz.eigen_vector[i].Y_SIZE = order;
		ritz.eigen_vector[i].V = malloc(sizeof(double) * ritz.eigen_vector[i].Y_SIZE);

	}
*/

//Algorithm run -----------------------------------------------

	
	gettimeofday(&timer_start, NULL);

	if (EXPLICIT_RESTART)
	{
		
		while(start_residu > TOLERENCE && iterations < CONVERGENCE_ITERATIONS)
		{	
			modified.V.X_SIZE = K_ITER;
			modified.V.Y_SIZE = M.X_SIZE;
			modified.H.X_SIZE = M.X_SIZE + 1;
			modified.H.Y_SIZE = M.X_SIZE + 1;

			modified.V.A = malloc(sizeof(double*)*modified.V.Y_SIZE);
			for (int i = 0; i < modified.V.Y_SIZE; ++i)
			{
				modified.V.A[i] = calloc(modified.V.X_SIZE, sizeof(double));
			}
			modified.H.A = malloc(sizeof(double*)*modified.H.Y_SIZE);
			for (int i = 0; i < modified.H.Y_SIZE; ++i)
			{
				modified.H.A[i] = calloc(modified.H.X_SIZE, sizeof(double));
			}

					
			Arnoldi_modified(modified,M,V,K_ITER,M.X_SIZE);
			modified.H.X_SIZE = M.X_SIZE;
			modified.H.Y_SIZE = M.X_SIZE;

		
			ritz_computation_by_arnoldi(modified.H, wr, wi, vr);
			
		
			sort_eigen(wr, vr, order);
			if (DEBUG_AFF)
			{
				printf("Matrix H: \n");
				Matrix_print(modified.H);

				printf("Vector V :\n");
				Matrix_print(modified.V);
				printf("Eigen values : \n");
		
				for (int i = 0; i < order; ++i)
				{
					printf("wr[%d] =  %f ::: ",i,wr[i]);
					for (int j = 0; j < order; ++j)
					{
						printf(" %f",vr[i*order +j] );
					}
					printf("\n");
			
				}
			
			}
		
		   	start_residu = residu_ritz(M, wr, vr, order);
		   	residu.V[iterations] = start_residu;


			
		   	if (iterations > 1)
		   	{
		   		
		   		convergence[iterations-1] = residu.V[iterations] - residu.V[iterations-1];
		   		convergence[iterations-1] = convergence[iterations-1] > 0 ? convergence[iterations-1] : -convergence[iterations-1];
		   		start_residu = convergence[iterations-1];

		   		printf("Convergence %d : %f\n",iterations, convergence[iterations-1] );
		   	}
	
	
		   	
		   	recompute_initial_vector_explicit(V, order, vr);
		   
		   	iterations ++;

		   	if (start_residu > TOLERENCE && iterations < CONVERGENCE_ITERATIONS)
		   	{
				for (int i = 0; i < modified.V.Y_SIZE; ++i)
				{
					free(modified.V.A[i]);
				}
			   	free(modified.V.A);
				for (int i = 0; i < modified.H.Y_SIZE; ++i)
				{
					free(modified.H.A[i]);
				}
				free(modified.H.A);
		   	}



		}
	}
	gettimeofday(&timer_stop, NULL);

		//Print results -------------------------------------------
  	if (ALL_PRINT)
  	{	
		printf("\n---------- EXIT ---------\n");
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
  	if (PERFORMANCE_MEASURE)
		{
			printf("Iterations, Time : %d %ld.%ld \n", iterations, timer_stop.tv_sec - timer_start.tv_sec, timer_stop.tv_usec - timer_start.tv_usec);
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