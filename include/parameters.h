#ifndef PARAMETERS_H
#define PARAMETERS_H 

#define RANDOM_OR_VERIFIED_3X3 1		// To generate random matrix and vector 
											// Random 1
											// 3x3 	  0 Predefined not working for MERAM
										// If aleatory generated matrix 
#define Matrix_size_X 5         			// Raw size
#define Matrix_size_Y 5					// Columns size 
#define MAX_VALUE_MATRIX 10				// Max Value in the matrix 
#define K_ITER 5						// N
#define DESIRED_EIGEN 5 					// Number of desired eigen
#define INITIAL_MATRIX_AFF 1			// Print the initial matrix    			  1 yes, 0 no 
#define INITIAL_VECTOR_AFF 1			// Print the initial vector    			  1 yes, 0 no 
#define H_MATRIX_AFF 1					// Print the result matrix H  			  1 yes, 0 no 
#define V_VECTOR_AFF 1					// Print the result vector V   			  1 yes, 0 no 
#define ALL_PRINT 1						// Print all matrix vectors 			  1 yes, 0 no
#define DEBUG_AFF 0 					// Dev check prints 					  1 yes, 0 no
#define PERFORMANCE_MEASURE	1			// Print the performance measures   	  1 yes, 0 no 
#define CONVERGENCE_PRINT 0 			// Print Convergence Value 				  1 yes, 0 no

#define ORTHOGONALITY_TEST 1			// Test if result vectors are orthogonal  1 yes, 0 no 
#define TOLERENCE 0.001					// Tolerence of residu 
#define CONVERGENCE_ITERATIONS 100		// Iterations to execute all the algorithm

#define EXPLICIT_RESTART 1 				// Run explicit restart method
#define MERAM 0							// Run MERAM
#define THREADS_NUMBER 2				// Threads number

#endif