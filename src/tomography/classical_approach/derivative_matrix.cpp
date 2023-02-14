# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "utils.hpp"


Utils::sparseMatrix Utils::getDerivativeMatrix(int n, int degree)
{
	sparseMatrix L;
    
    int elements = degree + 1;
		
	L.m = n;
    L.n = n - degree;
    L.nnz = elements * L.n;

    L.init();

	if (degree == 0)
	{
		for (int index = 0; index < L.nnz; index++)
		{
			L.i[index] = index;
			L.j[index] = index;
			L.v[index] = 1.0f;
		}
		
		return L;
	} 

	int * df = new int[elements]();	
	int * df1 = new int[elements + 1]();
	int * df2 = new int[elements + 1]();
	
	df[0] = -1; df[1] = 1;
	
	for (int index = 1; index < degree; index++)
	{
		for (int k = 0; k < elements; k++)
		{
			df2[k] = df[k];
			df1[k + 1] = df[k];

			df[k] = df1[k] - df2[k]; 
		}		 
	}
	
	for (int index = 0; index < L.n; index++)
	{
		for (int k = 0; k < elements; k++)
		{
			L.i[elements*index + k] = index;	
			L.j[elements*index + k] = index + k;
			L.v[elements*index + k] = df[k];
		}	
	}
	
	return L;
}


