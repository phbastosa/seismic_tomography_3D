# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <iostream>
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

    std::cout<<"Binary file " + path + " was successfully written."<<std::endl;

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

void Utils::sparse_lscg(int * iA, int * jA, float * vA, int n, int m, int nnz, float * B, float * x, int maxIt, float cgTol)
{    
    float a, b, qTq, rTr, rd;

    float * s = new float[n]();
    float * q = new float[n]();
    float * r = new float[m]();
    float * p = new float[m]();

    // s = d - A * x
    for (int row = 0; row < n; row++) 
        s[row] = B[row]; 

    // r = A' * s    
    for (int index = 0; index < nnz; index++) 
        r[jA[index]] += vA[index] * s[iA[index]];        

    // x = 0 & p = r 
    for (int col = 0; col < m; col++) 
    {
        x[col] = 0.0f;
        p[col] = r[col]; 
    }

    // q = A * p
    for (int index = 0; index < nnz; index++) 
        q[iA[index]] += vA[index] * p[jA[index]];        

    // # pragma acc enter data copyin(s[0:n], q[0:n], B[0:n])
    // # pragma acc enter data copyin(r[0:m], p[0:m], x[0:m])
    // # pragma acc enter data copyin(iA[0:nnz],jA[0:nnz],vA[0:nnz])

    for (int iteration = 0; iteration < maxIt; iteration++)
    {
        qTq = 0.0f;
        // # pragma acc parallel loop reduction(+:qTq) present(q[0:n])
        for (int row = 0; row < n; row++)                // q inner product
            qTq += q[row] * q[row];                      // qTq = q' * q

        rTr = 0.0f;
        // # pragma acc parallel loop reduction(+:rTr) present(r[0:m])
        for (int col = 0; col < m; col++)                // r inner product
            rTr += r[col] * r[col];                      // rTr = r' * r 

        a = rTr / qTq;                                   // a = (r' * r) / (q' * q)                    

        // # pragma acc parallel loop present(x[0:m], p[0:m])
        for (int col = 0; col < m; col++)                // model atualization
            x[col] += a * p[col];                        // x = x + a * p

        // # pragma acc parallel loop present(s[0:n], q[0:n])
        for (int row = 0; row < n; row++)                // s atualization  
            s[row] -= a * q[row];                        // s = s - a * q 

        rd = 0.0f;
        // # pragma acc parallel loop reduction(+:rd) present(r[0:m])
        for (int col = 0; col < m; col++)                // r inner product for division 
            rd += r[col] * r[col];                       // rd = r' * r

        // # pragma acc parallel loop present(r[0:m])
        for (int col = 0; col < m; col++)                // Zeroing r 
            r[col] = 0.0f;                               // r = 0, for multiplication

        // # pragma acc parallel loop present(iA[0:nnz],jA[0:nnz],vA[0:nnz],r[0:m],s[0:n])    
        for (int index = 0; index < nnz; index++)        // r atualization 
            r[jA[index]] += vA[index] * s[iA[index]];    // r = G' * s    
        
        rTr = 0.0f;    
        // # pragma acc parallel loop reduction(+:rTr) present(r[0:m])        
        for (int col = 0; col < m; col++)                // r inner product
            rTr += r[col] * r[col];                      // rTr = r' * r

        if (sqrtf(rd) < cgTol) break;                    // Convergence condition
        
        b = rTr / rd;                                    // b = (r' * r) / rd

        // # pragma acc parallel loop present(p[0:m],r[0:m]) 
        for (int col = 0; col < m; col++)          
            p[col] = r[col] + b * p[col];                // p = r + b * p 

        // # pragma acc parallel loop present(q[0:n])
        for (int row = 0; row < n; row++) 
            q[row] = 0.0f;                               // q = 0, for multiplication

        // # pragma acc parallel loop present(iA[0:nnz],jA[0:nnz],vA[0:nnz],q[0:n],p[0:m])
        for (int index = 0; index < nnz; index++) 
            q[iA[index]] += vA[index] * p[jA[index]];    // q = G * p           
    }
    
    // # pragma acc exit data delete(r[0:m], p[0:m])
    // # pragma acc exit data delete(s[0:n], q[0:n], B[0:n])
    // # pragma acc exit data delete(iA[0:nnz],jA[0:nnz],vA[0:nnz])
    // # pragma acc exit data copyout(x[0:m])

    delete[] s; delete[] q; delete[] r; delete[] p; 
}

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

float * Utils::movingAverageSmoothing(float * volume, int nx, int ny, int nz, int samples)
{
    int nPoints = nx*ny*nz;

    int init = (int) (samples / 2);
    
    float * smoothed = new float[nPoints]();

    for (int index = 0; index < nPoints; index++) 
        smoothed[index] = volume[index];

    for (int i = init; i < nz - init; i++)
    {
        for (int j = init; j < nx - init; j++)
        {
            for (int k = init; k < ny - init; k++)
            {
                float xs = 0.0f; 
                float ys = 0.0f;
                float zs = 0.0f;
                
                for (int s = 0; s < samples; s++)
                {
                    int p = s - init;

                    xs += volume[i + (j + p)*nz + k*nx*nz];
                    ys += volume[(i + p) + j*nz + k*nx*nz];
                    zs += volume[i + j*nz + (k + p)*nx*nz];
                }        

                xs *= 1.0f / samples;
                ys *= 1.0f / samples;
                zs *= 1.0f / samples;

                smoothed[i + j*nz + k*nx*nz] = (xs + ys + zs) / 3.0f;
            }
        }   
    }

    return smoothed;
}

float * Utils::gaussianFilterSmoothing(float * volume, int nx, int ny, int nz, float stdv, int samples)
{
    int init = samples / 2;
    int nPoints = nx * ny * nz;
    int nKernel = samples * samples * samples;

    float pi = 4.0f * atanf(1.0f); 

    float * kernel = new float[nKernel]();
    float * smoothed = new float[nPoints]();

    for (int i = 0; i < nPoints; i++) 
        smoothed[i] = volume[i];

    if (stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int y = -init; y <= init; y++)
        {
            for (int x = -init; x <= init; x++)
            {
                for (int z = -init; z <= init; z++)
                {          
                    int index = (z+init) + (x+init)*samples + (y+init)*samples*samples; 

                    kernel[index] = 1.0f / (2.0f*pi*stdv*stdv) * exp(-((x*x + y*y + z*z) / (2.0f*stdv*stdv)));
        
                    sum += kernel[index]; 
                }
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
    else 
    {
        int mid = samples / 2;

        kernel[mid + mid*samples + mid*samples*samples] = 1.0f;
    }

    for (int k = init; k < ny - init; k++)
    {   
        for (int j = init; j < nx - init; j++)
        {
            for (int i = init; i < nz - init; i++)
            {       
                float accum = 0.0f;
                
                for (int yk = 0; yk < samples; yk++)
                {      
                    for (int xk = 0; xk < samples; xk++)
                    {      
                        for (int zk = 0; zk < samples; zk++)
                        {   
                            int index = zk + xk*samples + yk*samples*samples;   
                            int partial = (i-init+zk) + (j-init+xk)*nz + (k-init+yk)*nx*nz; 

                            accum += volume[partial] * kernel[index];
                        }        
                    }
                }
                
                smoothed[i + j*nz + k*nx*nz] = accum;
            }
        }   
    }

    return smoothed;
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
