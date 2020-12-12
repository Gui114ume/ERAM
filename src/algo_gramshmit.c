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

	double ** h;
	double ** vout;
	double *vt;
	double *q;
	double *v_col;
	h = malloc(sizeof(double*)*m);
	vout = malloc(sizeof(double*)*m);
	vt = malloc(sizeof(double)*m);
	q = malloc(sizeof(double)*m);
	v_col = malloc(sizeof(double)*m);

	for (int i = 0; i < m; ++i)
	{
		h[i] = malloc(sizeof(double)*m);
		vout[i] = malloc(sizeof(double)*n);
		for (int j = 0; j < m; ++j)
		{
			h[i][j] = 0.;
		}
		for (int ncol = 0; ncol < n; ++ncol)
		{
			vout[i][ncol] = 0.;
			vt[ncol] = 0.;
			q[ncol] = 0.;
			v_col[ncol] = 0.;
		}
	}

	
	res.V.X_SIZE = n;
	res.V.Y_SIZE = m;
	
	
	res.H.X_SIZE = m;
	res.H.Y_SIZE = m;
	
	// INITIALIZATION -----------------------------------------
	
	 double module = norm_vector(v);
	
	for (int i = 0; i < res.V.Y_SIZE; ++i)
	{
		q[i] = v.V[i]/module;
	}
	for (int i = 0; i < res.V.Y_SIZE; ++i)
	{
		vout[i][0] = q[i];
	}



	// ALGORITHM ----------------------------------------------
	for (int j = 0; j < m-1; ++j)
	{
		
		dot_prod_mat_vec(vt, A, q);
		
		for (int i = 0; i < j+1; ++i)
		{
			for (int index = 0; index < res.V.Y_SIZE; ++index)
			{
				v_col[index] = vout[index][i];
			}

			h[i][j] = scalaire_vec_vec(vt, v_col, n); 

			for (int q = 0; q < n; ++q)
			{
				vt[q] = vt[q] - v_col[q] * h[i][j];
			}
		}
		h[j+1][j] = norme_vecteur_double(vt, n);

		if (h[j+1][j] > 0.0000000000001)
		{
		
			for (int i = 0; i < res.V.Y_SIZE; ++i)
			{
				q[i] = vt[i] / h[j+1][j];
				vout[i][j+1] = q[i];
			}
		}
		else
			break;
					
	}

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			res.H.A[i][j] = h[i][j];
		}
		for (int j = 0; j < n; ++j)
		{
			res.V.A[i][j] = vout[i][j];
		}
	}
	

	for (int i = 0; i < m; ++i)
	{
		free(h[i]);
		free(vout[i]);	
	}
	free(h);
	free(vout);
	free(vt);
	free(q);
	free(v_col);
	

	
}
