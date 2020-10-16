#include "algo_eram.h"

int find_max(double* module,
	     int nb_values)
{
  double max = module[0];
  int index = 0;
  for(int i = 1 ; i < nb_values ; i++)
    {
      if(module[i] > max)
	{
	  index = i;
	  max = module[i];
	}
    }
  return index;
}

Eigenvalues n_biggest(Eigenvalues eigenvalues,
		      int nb)
{
  printf("return the nb biggest eigenvalues\n");
  if( nb > eigenvalues.nb_values)
    {
      nb = eigenvalues.nb_values;
    }
  Eigenvalues res;
  res.nb_values = nb;
  res.real = malloc(sizeof(double) * nb);
  res.imag = malloc(sizeof(double) * nb);
  
  double module[eigenvalues.nb_values];
  int index_list[nb];
  for(int i = 0 ; i < eigenvalues.nb_values ; i++)
    {
      module[i] = sqrt(eigenvalues.real[i]*eigenvalues.real[i] + eigenvalues.imag[i]*eigenvalues.imag[i]);
    }
  int i = 0;
  int index = 0;
  while(i < nb)
    {
      index = find_max(module, eigenvalues.nb_values);
      index_list[i] = index;
      module[index] = -1;
      res.real[i] = eigenvalues.real[index];
      res.imag[i] = eigenvalues.imag[index];
      i++;
    }
  return res;
}

Eigenvalues compute_eigenvalues(Matrix_1d hessenberg)
{
  printf("compute eigenvalues\n");
  Eigenvalues eigenvalues;
  return eigenvalues;
}

Eigenvectors compute_eigenvectors(Matrix_1d   hessenberg,
				  Eigenvalues eigenvalues)
{
  printf("compute eigenvectors\n");
  Eigenvectors eigenvectors;
  return eigenvectors;
}

Vector create_initial_vector(Eigenvectors eigenvectors)
{
  printf("create initial vector\n");
  Vector vector;
  return vector;
}

int condition(double rs)
{
  printf("condition\n");
  double something;
  return something;
}

/* 
// to complete last

void Eram(A, v, n, m)
{
  printf("Eram\n");
  // choose parameter m and initial vector v
  int m = 10;
  int nb = 3;
  double tol = 1;
  double rs  = 2;
  Vector v;
  
  // iterate while we're not satisfied
  while(condition(rs) > tol)
    {
      // it does what it does
      Arnoldi_res arnoldi_res = Arnoldi_modified(A, v, n, m);

      // change matrix representation according to lapack convention
      Matrix_1d hessenberg= Matrix_to_Matrix_1d(arnoldi_res.H);
      
      // compute all eigenvalues
      Eigenvalues eigenvalues  = compute_eigenvalues(hessenberg);

      // extract the nb biggest eigenvalues (module)
      eigenvalues = n_biggest(eigenvalues,
			      nb);

      // compute eigenvectors associated with the nb eigenvalues
      Eigenvectors eigenvectors = compute_eigenvectors(hessenberg,
						      eigenvalues);
      // create the new initial vector with all the eigenvectors
      v = create_initial_vector(eigenvectors);
    }

  
    }*/
