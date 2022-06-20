# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <algorithm>

# include "utils.hpp"

int Utils::imin(int v1, int v2) { return !(v1 > v2) ? v1 : v2; }

int Utils::imax(int v1, int v2) { return !(v1 < v2) ? v1 : v2; }

float Utils::min(float v1, float v2) { return !(v1 > v2) ? v1 : v2; }

float Utils::max(float v1, float v2) { return !(v1 < v2) ? v1 : v2; }

float Utils::max3(float v1, float v2, float v3) { return max(v1, max(v2, v3)); }

float Utils::min3(float v1, float v2, float v3) { return min(v1, min(v2, v3)); }

float Utils::min4(float v1, float v2, float v3, float v4) { return min(v1, min(v2, min(v3, v4))); }

float * Utils::readBinaryFloat(std::string path, int n)
{
    std::ifstream file(path, std::ios::in);
    
    float * array = new float[n];

    if (file.is_open()) 
    {    
        file.read((char *) array, n * sizeof(float));
    }
    else
    {
        throw std::invalid_argument("Error: file could not be opened!");
    }

    file.close();    

    return array;
}

void Utils::writeBinaryFloat(std::string path, float *array, int n)
{
    std::ofstream file(path, std::ios::out);
    
    if (file.is_open()) 
    {    
        file.write((char *) array, n * sizeof(float));
    }
    else
    {
        throw std::invalid_argument("Error: file could not be opened!");
    }

    file.close();
}

std::string Utils::catchParameter(std::string target, std::string file)
{
    char spaces = ' ';
    char comment = '#';

    std::string line;
    std::string variable;

    std::ifstream parameters(file);

    if (parameters.is_open())
    {
        while (getline(parameters, line))
        {           
            if ((line.front() != comment) && (line.front() != spaces))        
            {
                if (line.find(target) == 0)
                {
                    for (int i = line.find("=")+2; i < line.size(); i++)
                    {    
                        if (line[i] == '#') break;
                        variable += line[i];            
                    }

                    break;
                }
            }                 
        }
        parameters.close();
    }        

    // Quality control for file paths

    if (variable.find('"') == 0)
    {
        remove(variable.begin(), variable.end(), '"');
    }
    else if (variable.find("[") == 0)
    {
        remove(variable.begin(), variable.end(), '[');
        remove(variable.begin(), variable.end(), ']');
    }

    variable.erase(remove(variable.begin(), variable.end(), ' '), variable.end());

    return variable;
}

std::vector<std::string> Utils::split(std::string s, char delimiter)
{
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) 
        tokens.push_back(token);
   
    return tokens;
}

bool Utils::str2bool(std::string s)
{
    bool b;

    std::for_each(s.begin(), s.end(), [](char & c) {c = ::tolower(c);});
    std::istringstream(s) >> std::boolalpha >> b;

    return b;
}

float * Utils::sparse_lscg(sparseMatrix A, float * B, int maxIt, float cgTol)
{    
    float a, b, qTq, rTr, rd;

    float * s = new float[A.n]();
    float * q = new float[A.n]();
    float * r = new float[A.m]();
    float * p = new float[A.m]();
    float * x = new float[A.m]();    // Linear system solution

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int row = 0; row < A.n; row++) 
        s[row] = B[row]; 

    // r = G' * s    
    for (int ind = 0; ind < A.nnz; ind++) 
        r[A.j[ind]] += A.v[ind] * s[A.i[ind]];        

    // p = r
    for (int col = 0; col < A.m; col++) 
        p[col] = r[col]; 

    // q = G * p
    for (int ind = 0; ind < A.nnz; ind++) 
        q[A.i[ind]] += A.v[ind] * p[A.j[ind]];        

    for (int iteration = 0; iteration < maxIt; iteration++)
    {
        qTq = 0.0f;
        for (int row = 0; row < A.n; row++)          // q inner product
            qTq += q[row] * q[row];                  // qTq = q' * q

        rTr = 0.0f;
        for (int col = 0; col < A.m; col++)          // r inner product
            rTr += r[col] * r[col];                  // rTr = r' * r 

        a = rTr / qTq;                               // a = (r' * r) / (q' * q)                    

        for (int col = 0; col < A.m; col++)          // model atualization
            x[col] += a * p[col];                    // x = x + a * p

        for (int row = 0; row < A.n; row++)          // s atualization  
            s[row] -= a * q[row];                    // s = s - a * q 

        rd = 0.0f;
        for (int col = 0; col < A.m; col++)          // r inner product for division 
            rd += r[col] * r[col];                   // rd = r' * r

        for (int col = 0; col < A.m; col++)          // Zeroing r 
            r[col] = 0.0f;                           // r = 0, for multiplication
        
        for (int ind = 0; ind < A.nnz; ind++)        // r atualization 
            r[A.j[ind]] += A.v[ind] * s[A.i[ind]];   // r = G' * s    

        rTr = 0.0f;                
        for (int col = 0; col < A.m; col++)          // r inner product
            rTr += r[col] * r[col];                  // rTr = r' * r

        if (sqrtf(rd) < cgTol) break;                // Convergence condition
        
        b = rTr / rd;                                // b = (r' * r) / rd

        for (int col = 0; col < A.m; col++)          
            p[col] = r[col] + b * p[col];            // p = r + b * p 

        for (int row = 0; row < A.n; row++) 
            q[row] = 0.0f;                         // q = 0, for multiplication

        for (int ind = 0; ind < A.nnz; ind++) 
            q[A.i[ind]] += A.v[ind] * p[A.j[ind]]; // q = G * p   
    }
    
    delete[] s; delete[] q; delete[] r; delete[] p;
    
    return x;
}

Utils::sparseMatrix Utils::firstOrderMatrixOperator(int order)
{
    sparseMatrix A;

    A.n = order;
    A.m = order;

    A.nnz = 3*(order-2) + 4;

    A.i = new int[A.nnz];
    A.j = new int[A.nnz];
    A.v = new float[A.nnz];
    
    for (int index = 0; index < order; index++)
    {
        if (index == 0)
        {
            A.i[index] = index;
            A.j[index] = index;
            A.v[index] = 1.0f;
            
            A.i[index + 1] = index;
            A.j[index + 1] = index + 1;
            A.v[index + 1] = -1.0f;
        }	
        else if ((index > 0) && (index < order-1))
        {
            A.i[2 + 3*(index - 1)] = index;
            A.j[2 + 3*(index - 1)] = index - 1;
            A.v[2 + 3*(index - 1)] = -1.0f;		

            A.i[3 + 3*(index - 1)] = index;
            A.j[3 + 3*(index - 1)] = index;
            A.v[3 + 3*(index - 1)] = 2.0f;		

            A.i[4 + 3*(index - 1)] = index;
            A.j[4 + 3*(index - 1)] = index + 1;
            A.v[4 + 3*(index - 1)] = -1.0f;
        }
        else	
        {
            A.i[A.nnz - 2] = index;
            A.j[A.nnz - 2] = index - 1;
            A.v[A.nnz - 2] = -1.0f;
            
            A.i[A.nnz - 1] = index;
            A.j[A.nnz - 1] = index;
            A.v[A.nnz - 1] = 1.0f;		
        }
    }

    return A;
}

Utils::sparseMatrix Utils::secondOrderMatrixOperator(int order)
{
	sparseMatrix A;
    
    A.n = order;
    A.m = order;

	A.nnz = 5*(order-4) + 14;

	A.i = new int[A.nnz];
	A.j = new int[A.nnz];
	A.v = new float[A.nnz];
	
	for (int index = 0; index < order; index++)
	{
		if (index == 0)
		{
			A.i[index] = index;
			A.j[index] = index;
			A.v[index] = 1.0f;
			
			A.i[index + 1] = index;
			A.j[index + 1] = index + 1;
			A.v[index + 1] = -2.0f;
			
			A.i[index + 2] = index;
			A.j[index + 2] = index + 2;
			A.v[index + 2] = 1.0f;			
		}
		else if (index == 1)
		{
			A.i[index + 2] = index;
			A.j[index + 2] = index - 1;
			A.v[index + 2] = -2.0f;
			
			A.i[index + 3] = index;
			A.j[index + 3] = index;
			A.v[index + 3] = 5.0f;
			
			A.i[index + 4] = index;
			A.j[index + 4] = index + 1;
			A.v[index + 4] = -4.0f;					

			A.i[index + 5] = index;
			A.j[index + 5] = index + 2;
			A.v[index + 5] = 1.0f;					
		}	
		else if ((index > 1) && (index < order-2))
		{
			A.i[2 + 5*(index - 1)] = index;
			A.j[2 + 5*(index - 1)] = index - 2;
			A.v[2 + 5*(index - 1)] = 1.0f;		

			A.i[3 + 5*(index - 1)] = index;
			A.j[3 + 5*(index - 1)] = index - 1;
			A.v[3 + 5*(index - 1)] = -4.0f;		

			A.i[4 + 5*(index - 1)] = index;
			A.j[4 + 5*(index - 1)] = index;
			A.v[4 + 5*(index - 1)] = 6.0f;		

			A.i[5 + 5*(index - 1)] = index;
			A.j[5 + 5*(index - 1)] = index + 1;
			A.v[5 + 5*(index - 1)] = -4.0f;

			A.i[6 + 5*(index - 1)] = index;
			A.j[6 + 5*(index - 1)] = index + 2;
			A.v[6 + 5*(index - 1)] = 1.0f;			
		}
		else if (index == order-2)	
		{
			A.i[A.nnz - 7] = index;
			A.j[A.nnz - 7] = index - 2;
			A.v[A.nnz - 7] = 1.0f;
			
			A.i[A.nnz - 6] = index;
			A.j[A.nnz - 6] = index - 1;
			A.v[A.nnz - 6] = -4.0f;
			
			A.i[A.nnz - 5] = index;
			A.j[A.nnz - 5] = index;
			A.v[A.nnz - 5] = 5.0f;					

			A.i[A.nnz - 4] = index;
			A.j[A.nnz - 4] = index + 1;
			A.v[A.nnz - 4] = -2.0f;					
		}
		else
		{
			A.i[A.nnz - 3] = index;
			A.j[A.nnz - 3] = index - 2;
			A.v[A.nnz - 3] = 1.0f;
			
			A.i[A.nnz - 2] = index;
			A.j[A.nnz - 2] = index - 1;
			A.v[A.nnz - 2] = -2.0f;					

			A.i[A.nnz - 1] = index;
			A.j[A.nnz - 1] = index;
			A.v[A.nnz - 1] = 1.0f;					
		}
	}

	return A;
}

float Utils::triLinearInterpolation(float c000, float c001, float c100, float c101, float c010, float c011, float c110, float c111, 
                                    float x0, float x1, float y0, float y1, float z0, float z1, float x, float y, float z)
{
    float xd = (x - x0) / (x1 - x0);
    float yd = (y - y0) / (y1 - y0);
    float zd = (z - z0) / (z1 - z0);

    float c00 = c000*(1 - xd) + c100*xd;    
    float c01 = c001*(1 - xd) + c101*xd;    
    float c10 = c010*(1 - xd) + c110*xd;    
    float c11 = c011*(1 - xd) + c111*xd;    

    float c0 = c00*(1 - yd) + c10*yd;
    float c1 = c01*(1 - yd) + c11*yd;

    return (c0*(1 - zd) + c1*zd);
}
