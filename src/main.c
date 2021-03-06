#include "types.h"
#include "algo_gramshmit.h"
#include "parameters.h"
#include "ritz_computation_by_arnoldi.h"
#include <pthread.h>

pthread_barrier_t   barrier;


void * ERAM_thread(void * args)
{
	Threads_str *str  = (Threads_str*)args;
	Arnoldi_res modified;
	double * wr;
	double * wi;
	double * vr;
	int order = str->M.X_SIZE;

	wr = malloc(sizeof(double)* order);
	wi = malloc(sizeof(double)* order);
	vr = malloc(sizeof(double)* order * order);

	modified.V.X_SIZE = K_ITER;
	modified.V.Y_SIZE = Matrix_size_Y;
	modified.H.X_SIZE = Matrix_size_X + 1;
	modified.H.Y_SIZE = Matrix_size_Y + 1;

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
	
	if (str->vector_updated == 0)
	{
		for (int i = 0; i < str->V.Y_SIZE; ++i)
		{
			str->V.V[i] = rand() % MAX_VALUE_MATRIX;
		}
		str->vector_updated++;
	}



	Arnoldi_modified(modified,str->M,str->V,K_ITER,str->M.X_SIZE);
	modified.H.X_SIZE = str->M.X_SIZE;
	modified.H.Y_SIZE = str->M.X_SIZE;

	ritz_computation_by_arnoldi(modified.H, wr, wi, vr);
	

	sort_eigen(wr, vr, order, DESIRED_EIGEN);
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

   	str->residu = residu_ritz(str->M, wr, vr, order);
   	recompute_initial_vector_explicit(str->V, order, vr);
 
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

	free(wr);
	free(wi);
	free(vr);


	pthread_exit(NULL);		
}

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
	Threads_str threads_args[THREADS_NUMBER];

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
			V.V[i] = rand() % MAX_VALUE_MATRIX;	
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
			
		
			sort_eigen(wr, vr, order, DESIRED_EIGEN);
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
		   		
		   		convergence[iterations-1] = 1 - (residu.V[iterations] / residu.V[iterations-1]);
		   		convergence[iterations-1] = convergence[iterations-1] > 0 ? convergence[iterations-1] : -convergence[iterations-1];
		   		start_residu = convergence[iterations-1];
		   		if (CONVERGENCE_PRINT)
		   		{
		   			printf("Convergence %d : %f\n",iterations, convergence[iterations-1] );
		   		}
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
	if (MERAM)
	{
		for (int id = 0; id < THREADS_NUMBER; ++id)
		{
		threads_args[id].vector_updated = 0;
		threads_args[id].V.Y_SIZE = Matrix_size_Y;
		threads_args[id].M.X_SIZE = M.X_SIZE;
		threads_args[id].M.Y_SIZE = M.Y_SIZE;
		threads_args[id].V.V = malloc(sizeof(double)*threads_args[id].V.Y_SIZE);
		threads_args[id].M.A = malloc(sizeof(double*)*threads_args[id].M.Y_SIZE);
		for (int i = 0; i < threads_args[id].M.Y_SIZE; ++i)
		{
			threads_args[id].M.A[i] = malloc(sizeof(double)*threads_args[id].M.X_SIZE);
			for (int j = 0; j < threads_args[id].M.X_SIZE; ++j)
			{
				threads_args[id].M.A[i][j] = M.A[i][j];
			}
		}

		}

		pthread_t tid[THREADS_NUMBER];
		int end = 0;
		double convergence_thread[THREADS_NUMBER][CONVERGENCE_ITERATIONS];
		double residu_thread[THREADS_NUMBER][CONVERGENCE_ITERATIONS];

		pthread_barrier_init (&barrier, NULL, THREADS_NUMBER);
		while(!end)
		{
			for (int i = 0; i < THREADS_NUMBER; ++i)
			{
				pthread_create(&tid[i], NULL, ERAM_thread, &threads_args[i] );
			}

			for (int i = 0; i < THREADS_NUMBER; ++i)
			{
				pthread_join(tid[i], NULL);
		
			}
			for (int i = 0; i < THREADS_NUMBER; ++i)
			{
				if (iterations > 1)
		   		{
		   		
		   			convergence_thread[i][iterations-1] = 1 - (residu_thread[i][iterations] / residu_thread[i][iterations-1]);
		   			convergence_thread[i][iterations-1] = convergence_thread[i][iterations-1] > 0 ? convergence_thread[i][iterations-1] : -convergence_thread[i][iterations-1];
		   			if (CONVERGENCE_PRINT)
		   			{
		   				printf("Convergence thread %d, %d : %f\n",i, iterations, convergence_thread[i][iterations-1] );
		   			}
		   			if (convergence_thread[i][iterations-1] < TOLERENCE && iterations > CONVERGENCE_ITERATIONS)
		   			{
		   				end = 1;
		   			}
		  	 	}
			}
			iterations++;

		}

		for (int i = 0; i < THREADS_NUMBER; ++i)
		{
			for (int j = 0; j < threads_args[i].M.Y_SIZE; ++j)
			{
				free(threads_args[i].M.A[j]);
			}
			free(threads_args[i].V.V);
			free(threads_args[i].M.A);
	
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
			long int sec = timer_stop.tv_sec - timer_start.tv_sec;
			long int usec = timer_stop.tv_usec - timer_start.tv_usec;
			if (sec)
			{
				usec = usec < 0 ? -usec : usec;
			}
			printf("Iterations, Time : %d %ld.%ld \n", iterations, sec, usec);
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