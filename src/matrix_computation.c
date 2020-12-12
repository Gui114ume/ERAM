#include "matrix_computation.h"	

void Matrix_print (Matrix M)
{
	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < M.X_SIZE; ++j)
		{
			printf("%.2f ",M.A[i][j] );
		}
		printf("\n");
	}
	printf("\n");

}
void Vector_print (Vector V)
{
	for (int i = 0; i < V.Y_SIZE; ++i)
	{
		
		printf("%.2f \n",V.V[i] );
		
	}
	printf("\n");

}

Matrix dot_product_matrix(Matrix A, Matrix B)
{
	Matrix res;
	res.X_SIZE = B.X_SIZE;
	res.Y_SIZE = A.Y_SIZE;

	res.A = malloc(sizeof(double*)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		res.A[i] = malloc(sizeof(double)*res.X_SIZE);
	}

	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		for (int j = 0; j < res.X_SIZE; ++j)
		{
			for (int m = 0; m < A.X_SIZE; m++)
			{	
				res.A[i][j] += A.A[i][m] * B.A[m][j];
				//printf("%d %d += %d %d * %d %d \n",i,j,i,m,m,j );

			}
		}
	}


	return res;
}

void dot_prod_matrix_vector_from_matrix(Vector vout, Matrix A, Matrix B, int column)
{	
	if (A.X_SIZE != B.Y_SIZE)
	{
		printf("Dot prod matrix vector not possible %d != %d \n",A.X_SIZE, B.Y_SIZE );
	}
	
	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		vout.V[i] = 0;
		for (int m = 0; m < A.X_SIZE; ++m)
		{
			vout.V[i] += A.A[i][m] * B.A[m][column];
		}
	}

}

double dot_product_vector(Vector A, Vector B)
{
	if (A.Y_SIZE != B.Y_SIZE)
	{
		printf("Vector prod error. Check Vectors lenght\n");
		return 0;
	}
	double res = 0.;
	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		res += A.V[i] * B.V[i];
	}

	return res;
}
double dot_product_vector_matrix_column(Vector A, Matrix B, int column)
{
	if (A.Y_SIZE != B.Y_SIZE)
	{
		printf("Vector prod error. Check Vectors lenght\n");
		printf("Ay : %d By : %d\n",A.Y_SIZE,B.Y_SIZE);
		return 0;
	}
	double res = 0.;
	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		res += A.V[i] * B.A[i][column];
	}

	return res;
}
Vector prod_matrix_column_coef(Matrix A, int column, double k)
{
	Vector res;
	res.Y_SIZE = A.Y_SIZE;
	res.V = malloc(sizeof(double)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		res.V[i] = A.A[i][column] * k;
	}

	return res;
}
void vector_minus_vector(Vector A, Vector B)
{

	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		A.V[i] = A.V[i] - B.V[i];
	}
}

double norm_vector(Vector A)
{	
	double dot_prod = dot_product_vector(A,A);
	double norm = sqrt(dot_prod);
	//for (int i = 0; i < A.Y_SIZE; ++i)
	//{
	//	A.V[i] /= norm;
	//}
	return norm;
}
void normalize_vector(Vector A)
{	
	double dot_prod = dot_product_vector(A,A);
	double norm = sqrt(dot_prod);
	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		A.V[i] /= norm;
	}
}

int check_orthogolality(Matrix A)  // Each column is a vector
{
	int res  = 0;
	for (int i = 0; i < A.X_SIZE; ++i)
	{
		for (int j = i+1; j < A.X_SIZE; ++j)
		{
			for (int col = 0; col < A.Y_SIZE; ++col)
			{
				res += A.A[i][col] * A.A[j][col];
			}
		}
		if (res)
		{
			return -1;
		}
	}


	return res;
}

void from2d_to_1d(double * M_oneD, Matrix M)
{	
	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < M.X_SIZE; ++j)
		{
			M_oneD[i*M.X_SIZE+j] = M.A[i][j];
		}
	}
}

Matrix coeff_dot_identy(double k, int order)
{
	Matrix m;
	m.X_SIZE = order;
	m.Y_SIZE = order;
	m.A = malloc(sizeof(double*)*m.Y_SIZE);
	for (int i = 0; i < m.Y_SIZE; ++i)
	{
		m.A[i] = malloc(sizeof(double)*m.X_SIZE);
	}
	
	for (int i = 0; i < m.X_SIZE; ++i)
	{
		m.A[i][i] = k;
	}

	return m;
}

Matrix matrix_minus_matrix(Matrix A, Matrix B)
{	
	if (A.X_SIZE != B.X_SIZE || A.Y_SIZE != B.Y_SIZE)
	{
		printf("Error matrix minus matrix \n");
		exit(1);
	}
	Matrix m;
	m.X_SIZE = A.X_SIZE;
	m.Y_SIZE = A.Y_SIZE;
	m.A = malloc(sizeof(double*)*m.Y_SIZE);
	for (int i = 0; i < m.Y_SIZE; ++i)
	{
		m.A[i] = malloc(sizeof(double)*m.X_SIZE);
	}
	
	for (int i = 0; i < m.Y_SIZE; ++i)
	{
		for (int j = 0; j < m.X_SIZE; ++j)
		{	
			m.A[i][j] = A.A[i][j] - B.A[i][j];
		}
	}

	return m;

}

Vector dot_Matrix_Vector(Matrix A, double* vr, int order, int index)
{
	Vector res;
	res.Y_SIZE = A.Y_SIZE;
	res.V = malloc(sizeof(double)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		res.V[i] = 0;
		for (int m = 0; m < A.X_SIZE; ++m)
		{
			res.V[i] += A.A[i][m] * vr[index * order + m];
		}
	}

	return res;
}

double vector_norm(Vector v)
{
	double res = 0.;
	for (int i = 0; i < v.Y_SIZE; ++i)
	{
		res += v.V[i] * v.V[i];
	}
	res = sqrt(res);
	return res; 
}

double norm_matrix_column(Matrix A, int column)
{
	double res = 0.;
	for (int i = 0; i < A.Y_SIZE; ++i)
	{
		res+= A.A[i][column] * A.A[i][column];
	}
	res = sqrt(res);
	return res;
}

void sort_eigen(double* wr, double* vr, int order)
{
	
	int index[order];
	int sorted = 0;
	double tmp_eingen = 0.;
	int tmp_index; 
	double tmp_eingen_vec[order*order];

	
	for (int i = 0; i < order; ++i)
	{
		index[i] = i;
	}

	
	// Sorting eigenvalues and vectors ------------------------
	if (order > 2)
	{
		while(!sorted)
		{
			sorted = 1;
			for (int i = 0; i < order - 1; ++i)
			{
				if (abs(wr[i]) < abs(wr[i+1]))
				{
					tmp_eingen = wr[i];
					wr[i] = wr[i+1];
					wr[i+1] =tmp_eingen;
					sorted = 0;
					tmp_index = index[i];
					index[i] = index[i+1];
					index[i+1] = tmp_index;
				}
			}
		}
		
		
	}
	else
	{
		if (order == 1)
		{
			printf("Invalid matrix order \n");
			exit(1);
		}
		if (abs(wr[0] > abs(wr[1])))
		{
			tmp_eingen = wr[0];
			wr[0] = wr[1];
			wr[1] = tmp_eingen;

		}
	}
	
	for (int i = 0; i < order*order; ++i)
	{
		tmp_eingen_vec[i] = vr[i];
	}

	// Ordering vectors -------------
	for (int i = 0; i < order; ++i)
	{
		for (int j = 0; j < order; ++j)
		{
			vr[i*order + j] = tmp_eingen_vec[index[i]*order + j];
		}
	}
	
}

void recompute_initial_vector_explicit(Vector initial, int order, double *vr)
{
	if (order != initial.Y_SIZE)
	{
		printf("Unable to compute initial vector \n");
		exit(1);
	}
	double tmp[order];
	for (int i = 0; i < order; ++i)
	{
		tmp[i] = initial.V[i];
	}
	for (int i = 0; i < order; ++i)
	{
		initial.V[i] = 0;
		for (int j = 0; j < order; ++j)
		{
			initial.V[i] += tmp[i] * vr[order * i + j];
		}
	}
}

void QRfactorisation(Matrix M, Matrix Q, Matrix R)
{	
	Matrix tmp_a;
	Matrix tmp_b;
	double r = 0.;
	tmp_a.Y_SIZE = M.Y_SIZE;
	tmp_a.X_SIZE = M.X_SIZE;
	tmp_b.Y_SIZE = M.Y_SIZE;
	tmp_b.X_SIZE = M.X_SIZE;

	tmp_a.A = malloc(sizeof(double*) * tmp_a.Y_SIZE);
	tmp_b.A = malloc(sizeof(double*) * tmp_b.Y_SIZE);
	for (int i = 0; i < tmp_a.Y_SIZE; ++i)
	{
		tmp_a.A[i] = malloc(sizeof(double)*tmp_a.X_SIZE);
		tmp_b.A[i] = malloc(sizeof(double)*tmp_b.X_SIZE);
	}


	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < M.X_SIZE; ++j)
		{
			Q.A[i][j] = M.A[i][j];
		}
	}
	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < i; ++j)
		{	
			r = 0.;
			for (int raw = 0; raw < M.Y_SIZE; ++raw)
			{
				tmp_a.A[raw][0] = Q.A[raw][j];
				tmp_b.A[raw][0] = M.A[raw][i];
				r += tmp_a.A[raw][0] * tmp_b.A[raw][0];
			}

			R.A[j][i] = r;
			for (int i = 0; i < tmp_a.Y_SIZE; ++i)
			{
				tmp_a.A[i][0] *= r;
			}
			for (int raw = 0; raw < Q.Y_SIZE; ++raw)
			{
				Q.A[raw][i] -= tmp_a.A[raw][j];
			}
		}
		R.A[i][i] = norm_matrix_column(Q,i);
		for (int raw = 0; raw < Q.Y_SIZE; ++raw)
		{	
			if (Q.A[raw][i])
			{
				Q.A[raw][i] /= R.A[i][i];
			}
		}
	}
	
	for (int i = 0; i < tmp_a.Y_SIZE; ++i)
	{
		free(tmp_a.A[i]);
		free(tmp_b.A[i]);
	}
	free(tmp_a.A);
	free(tmp_b.A);

}

void recompute_initial_vector_implicit(Matrix H, Matrix V, Ritz_eigen ritz, int shift_p)
{
	Matrix Q,R,Q_transp;
	Q.X_SIZE = H.X_SIZE;
	Q.Y_SIZE = H.Y_SIZE;
	R.X_SIZE = H.X_SIZE;
	R.Y_SIZE = H.Y_SIZE;
	Q.A = malloc(sizeof(double*)*Q.Y_SIZE);
	R.A = malloc(sizeof(double*)*R.Y_SIZE);
	for (int i = 0; i < H.Y_SIZE; ++i)
	{
		Q.A[i] = calloc(Q.X_SIZE, sizeof(double));
		R.A[i] = calloc(R.X_SIZE, sizeof(double));
	}
	Matrix shifted_H;
	shifted_H.X_SIZE = H.X_SIZE;
	shifted_H.Y_SIZE = H.Y_SIZE;
	shifted_H.A = malloc(sizeof(double*)*shifted_H.Y_SIZE);
	for (int i = 0; i < shifted_H.Y_SIZE; ++i)
	{
		shifted_H.A[i] = calloc(shifted_H.X_SIZE, sizeof(double));
	}

	for (int s = 0; s < shift_p; ++s)
	{
		for (int i = 0; i < H.Y_SIZE; ++i)
		{
			for (int j = 0; j < H.X_SIZE; ++j)
			{
				shifted_H.A[i][j] = H.A[i][j];
				if (i == j && i <= shift_p)
				{
					shifted_H.A[i][j] -= ritz.eigen_value[s];
				}
			}
		}
		QRfactorisation(shifted_H, Q, R);
		
		Q_transp = transposed(Q);
	
		H = dot_product_matrix(Q_transp, H);
	
		H = dot_product_matrix(H,Q);


		V = dot_product_matrix(V,Q);

	}
	

	for (int i = 0; i < H.Y_SIZE; ++i)
	{
		free(Q.A[i]);
		free(R.A[i]);
	}
	free(Q.A);
	free(R.A);
}

void vector_dot_coeff(Vector V, double k)
{
	for (int i = 0; i < V.Y_SIZE; ++i)
	{
		V.V[i] *= k;
	}
}


Matrix transposed(Matrix M)
{
	Matrix res;
	res.X_SIZE = M.Y_SIZE;
	res.Y_SIZE = M.X_SIZE;
	res.A = malloc(sizeof(double*)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		res.A[i] = calloc(res.X_SIZE, sizeof(double));
	}
	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < M.X_SIZE; ++j)
		{
			res.A[j][i] = M.A[i][j];
		}
	}
	return res;
}