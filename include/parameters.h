#ifndef PARAMETERS_H
#define PARAMETERS_H 

#define RANDOM_OR_VERIFIED_3X3 1	// To generate random matrix and vector 
										// Random 1
										// 3x3 	  0
									// If aleatory generated matrix 
#define Matrix_size_X 10         		// Raw size
#define Matrix_size_Y 10				// Columns size 
#define MAX_VALUE_MATRIX 10				// Max Value in the matrix 
#define K_ITER 10						// N

#define INITIAL_MATRIX_AFF 1		// Print the initial matrix    			  1 yes, 0 no 
#define INITIAL_VECTOR_AFF 1		// Print the initial vector    			  1 yes, 0 no 
								
#define H_MATRIX_AFF 1				// Print the result matrix H  			  1 yes, 0 no 
#define V_VECTOR_AFF 1				// Print the result vector V   			  1 yes, 0 no 
#define ALL_PRINT 1					// Print all matrix vectors 
#define PERFORMANCE_MEASURE	0		// Print the performance measures   	  1 yes, 0 no 
#define ORTHOGONALITY_TEST 1		// Test if result vectors are orthogonal  1 yes, 0 no 

#define TOLERENCE 0.00000000001		// Tolerence of residu 
#define CONVERGENCE_ITERATIONS 300	// Iterations to execute all the algorithm

#define EXPLICIT_RESTART 1 			// Run explicit restart method
#define IMPLICIT_RESTART 0 			// Run implicit restart method

#endif







