#include "ritz_computation_by_arnoldi.h"
#include "types.h"
/*
Parameters
[in]	JOBVL	
          JOBVL is CHARACTER*1
          = 'N': left eigenvectors of A are not computed;
          = 'V': left eigenvectors of A are computed.
[in]	JOBVR	
          JOBVR is CHARACTER*1
          = 'N': right eigenvectors of A are not computed;
          = 'V': right eigenvectors of A are computed.
[in]	N	
          N is INTEGER
          The order of the matrix A. N >= 0.
[in,out]	A	
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the N-by-N matrix A.
          On exit, A has been overwritten.
[in]	LDA	
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[out]	WR	
          WR is DOUBLE PRECISION array, dimension (N)
[out]	WI	
          WI is DOUBLE PRECISION array, dimension (N)
          WR and WI contain the real and imaginary parts,
          respectively, of the computed eigenvalues.  Complex
          conjugate pairs of eigenvalues appear consecutively
          with the eigenvalue having the positive imaginary part
          first.
[out]	VL	
          VL is DOUBLE PRECISION array, dimension (LDVL,N)
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues.
          If JOBVL = 'N', VL is not referenced.
          If the j-th eigenvalue is real, then u(j) = VL(:,j),
          the j-th column of VL.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
          u(j+1) = VL(:,j) - i*VL(:,j+1).
[in]	LDVL	
          LDVL is INTEGER
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N.
[out]	VR	
          VR is DOUBLE PRECISION array, dimension (LDVR,N)
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues.
          If JOBVR = 'N', VR is not referenced.
          If the j-th eigenvalue is real, then v(j) = VR(:,j),
          the j-th column of VR.
          If the j-th and (j+1)-st eigenvalues form a complex
          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
          v(j+1) = VR(:,j) - i*VR(:,j+1).
[in]	LDVR	
          LDVR is INTEGER
          The leading dimension of the array VR.  LDVR >= 1; if
          JOBVR = 'V', LDVR >= N.
[out]	WORK	
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]	LWORK	
          LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= max(1,3*N), and
          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
          performance, LWORK must generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]	INFO	
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors have been computed;
                elements i+1:N of WR and WI contain eigenvalues which
                have converged.
 */   
void eigenvalues_print(int order, double * wr, double * wi)
{
	/*printf("\nEigen values : \n");
	for (int i = 0; i < order; ++i)
	{
		printf("(%f, %f)\n",wr[i], wi[i]);
	}
	*/
	printf("\nEigen values real part : \n");
	for (int i = 0; i < order; ++i)
	{
		printf("%.2f \n",wr[i]);
	}
}
void eigenvector_print(char * vector_type, int order, double * wi, double * vector, int ldv)
{
	/*
	printf("\n%s : \n",vector_type);
   for( int i = 0; i < order; i++ ) {
     for (int j = 0; j < order; ++j)
     {
     	printf( " (%f,%f)", vector[i+j*ldv], vector[i+(j+1)*ldv] );
     }
      printf( "\n" );
   }
   */
	printf("\n%s real part : \n",vector_type);
   for( int i = 0; i < order; i++ ) {
     for (int j = 0; j < order; ++j)
     {
     	printf( " %.2f ", vector[i+j*ldv]);
     }
      printf( "\n" );
   }

}

            
void ritz_computation_by_arnoldi(Matrix M, double * wr, double * wi, double * vr)
{	
	int order = M.X_SIZE;
	int lda = order;	
	int ldvl = order;    											 // leading dimension of the array VL
	int ldvr = order;    	  										 // leading dimension of the array VR
	int lwork = order;      										 // dimension (MAX(1,LWORK))

	double * M_oneD = malloc(sizeof(double)*M.X_SIZE * M.Y_SIZE); 	 // matrice M en une dimension
									//wr     					 	 // partie reele 
									//wi      					 	 // partie imaginaire
	double vl[ldvl * order];     									 // left vector dimension (LDVL,N)
									//vr 		  	 				 // right vector dimension (LDVR,N)
	double *work;			 									 	 // the optimal LWORK dimension (MAX(1,LWORK))
	double wkopt;
	int info;		 												 // return status
	
	from2d_to_1d(M_oneD, M);
	//Allocate the optimal workspace --------------------------
	lwork = -1;
	dgeev("N","V", &order,M_oneD,
                &lda, wr, wi, vl, &ldvl,
                vr, &ldvr, &wkopt, &lwork, &info);
	lwork = (int)wkopt;
    work = (double*)calloc(sizeof(double),lwork);
    if (work == NULL)
    {
    	printf("Allocate work error\n");
    	exit(0);
    }


    //Compute eigen values and vectors -----------------------
    dgeev( "N", "V", &order, M_oneD, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info );

    if (info !=0 )
    {
    	printf("Error computing eigen values and vectors\n");
    	exit(1);
    }


	//free(work);
   	
   	free(M_oneD);

}
/*
double residu_ritz(Matrix M, double * wr, double * vr, int order)
{	
	double r[order];
	double r_tot = 0.;
	Matrix lambdaI;
	Matrix A_minus_lambdaI;
	Vector inside_norm;
	// Computing residu ---------------------------------------

	for (int i = 0; i < order; ++i)
	{	
		lambdaI = coeff_dot_identy(wr[i], order);
		A_minus_lambdaI = matrix_minus_matrix(M,lambdaI);

		inside_norm = dot_Matrix_Vector(A_minus_lambdaI, vr, order, i);
		r[i] = vector_norm(inside_norm);

		for (int f = 0; f <  lambdaI.Y_SIZE; ++f)
		{
			free(lambdaI.A[f]);
			free(A_minus_lambdaI.A[f]);
		}
			free(inside_norm.V);
			free(lambdaI.A);
			free(A_minus_lambdaI.A);
	
	}

	for (int i = 0; i < order; ++i)
	{
		r_tot += r[i] * r[i];
	}


	return sqrt(r_tot);

}
*/


double residu_ritz(Matrix M, double * wr, double * vr, int order)
{	
	double r[order];
	double r_tot = 0.;
	double ass_vector[order];
	Matrix A_minus_lambdaI;
	Matrix lambdaI;
	double ident_ass_vec[order];
	double mat_ass_vect[order];
	double inside_norm[order];

	// Computing residu ---------------------------------------

	for (int i = 0; i < order; ++i)
	{	
		lambdaI = coeff_dot_identy(wr[i], order);

		for (int index = 0; index < order; ++index)
		{
			ass_vector[index] = vr[i*order + index];
		}

		A_minus_lambdaI = matrix_minus_matrix(M,lambdaI);
		dot_prod_mat_vec(ident_ass_vec, A_minus_lambdaI, ass_vector);

		dot_prod_mat_vec(mat_ass_vect, M, ass_vector);
		for (int index = 0; index < order; ++index)
		{
			inside_norm[index] =  mat_ass_vect[index] - ident_ass_vec[index];
		}

		r[i] = norme_vecteur_double(inside_norm, order);

		for (int f = 0; f <  lambdaI.Y_SIZE; ++f)
		{	
			free(lambdaI.A[f]);
			free(A_minus_lambdaI.A[f]);
		}
			free(A_minus_lambdaI.A);
			free(lambdaI.A);
	
	}

	for (int i = 0; i < order; ++i)
	{
		r_tot += r[i] * r[i];
	}


	return sqrt(r_tot);

}