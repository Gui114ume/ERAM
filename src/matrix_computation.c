#include "matrix_computation.h"	

void Matrix_print (Matrix M)
{
	for (int i = 0; i < M.Y_SIZE; ++i)
	{
		for (int j = 0; j < M.X_SIZE; ++j)
		{
			printf("%f ",M.A[i][j] );
		}
		printf("\n");
	}
	printf("\n");

}

void Matrix_1d_print(Matrix_1d M)
{
  for (int i = 0; i < M.Y_SIZE; ++i)
    {
      for (int j = 0; j < M.X_SIZE; ++j)
	{
	  printf("%f ",M.A[i * M.X_SIZE + j] );
	}
      printf("\n");
    }
  printf("\n");
}

void Vector_print (Vector V)
{
	for (int i = 0; i < V.Y_SIZE; ++i)
	{
		
		printf("%f \n",V.V[i] );
		
	}
	printf("\n");

}

Matrix_1d Matrix_to_Matrix_1d(Matrix M)
{
  Matrix_1d res;
  res.X_SIZE = M.X_SIZE;
  res.Y_SIZE = M.Y_SIZE;
  res.A = malloc(sizeof(double) * res.X_SIZE * res.Y_SIZE);
  for(int i = 0 ; i < M.Y_SIZE ; ++i)
    {
      for(int j = 0 ; j < M.X_SIZE ; ++j)
	{
	  res.A[i * res.X_SIZE + j] = M.A[i][j];
	}
    }
  
  return res;
}

Matrix dot_product_matrix(Matrix A, Matrix B)
{
	Matrix res;
	res.X_SIZE = B.X_SIZE;
	res.Y_SIZE = A.Y_SIZE;
	if (A.X_SIZE != B.Y_SIZE)
	{
		printf("Dot product not defined : A.y != B.x\n");
		return res;	
	}

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
Vector dot_prod_matrix_vector_from_matrix(Matrix A, Matrix B, int column)
{	
	if (A.X_SIZE != B.Y_SIZE)
	{
		printf("Dot prod matrix vector not possible %d != %d \n",A.X_SIZE, B.Y_SIZE );
	}
	Vector res;
	res.Y_SIZE = A.Y_SIZE;
	res.V = malloc(sizeof(double)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		for (int m = 0; m < A.X_SIZE; ++m)
		{
			res.V[i] += A.A[i][m] * B.A[m][column];
		}
	}

	return res;

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
Vector vector_minus_vector(Vector A, Vector B)
{
	Vector res;
	res.Y_SIZE = A.Y_SIZE;
	res.V = malloc(sizeof(double)*res.Y_SIZE);
	for (int i = 0; i < res.Y_SIZE; ++i)
	{
		res.V[i] = A.V[i] - B.V[i];
	}
	return res;
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

