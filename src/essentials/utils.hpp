# ifndef UTILS_HPP
# define UTILS_HPP

# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <algorithm>

/* 
Struct to define a sparse matrix:
    n - number of rows
    m - number of cols
    nnz - number of non zero elements
    i - row index
    j - col index
    v - values 
*/ 
typedef struct     
{
    int * i;       // Rows indexes
    int * j;       // Cols indexes 
    float * v;     // Value

    int n;         // Rows number
    int m;         // Cols number
    int nnz;       // Non-zero elements

    /* Function to initialize the vectorial components of sparse matrix */
    void init(int nrows, int ncols, int nonZeros)
    {
        n = nrows;
        m = ncols;
        nnz = nonZeros;

        i = new int[nnz]();  
        j = new int[nnz](); 
        v = new float[nnz](); 
    }

    /* Function to delete the vectorial components of sparse matrix */
    void erase()
    {
        delete[] i;
        delete[] j;
        delete[] v;
    }

} sparseMatrix;

/* Reads and returns a binary float */
float * readBinaryFloat(std::string path, int n)
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

/* Writes a binary float */
void writeBinaryFloat(std::string path, float *array, int n)
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

/* Finds each parameter in text file */
std::string catchParameter(std::string target, std::string file)
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
            std::cout<<line<<std::endl;        
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

/* Function to calculate minimum value between two float inputs */
float min(float v1, float v2) 
{ 
    return !(v1 > v2) ? v1 : v2; 
}

/* Function to calculate maximum value between two float inputs */
float max(float v1, float v2) 
{ 
    return !(v1 < v2) ? v1 : v2; 
}

/* Function to calculate minimum value between two integer inputs */
int imin(int v1, int v2) 
{ 
    return !(v1 > v2) ? v1 : v2; 
}

/* Function to calculate maximum value between two integer inputs */
int imax(int v1, int v2) 
{ 
    return !(v1 < v2) ? v1 : v2; 
}
    
/* Function to calculate minimum value between three float inputs */
float min3(float v1, float v2, float v3) 
{ 
    return min(v1, min(v2, v3)); 
}

/* Function to calculate maximum value between three float inputs */
float max3(float v1, float v2, float v3) 
{ 
    return max(v1, max(v2, v3)); 
}

/* Function to calculate minimum value between four float inputs */
float min4(float v1, float v2, float v3, float v4) 
{ 
    return min(v1, min(v2, min(v3, v4))); 
}

/* Function to convert string to boolean */
bool str2bool(std::string s)
{
    bool b;

    std::for_each(s.begin(), s.end(), [](char & c) {c = ::tolower(c);});
    std::istringstream(s) >> std::boolalpha >> b;

    return b;
}

/* Function to separete values with a delimiter */
std::vector<std::string> split(std::string s, char delimiter)
{
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) 
        tokens.push_back(token);
   
    return tokens;
}

/* Function to compute a sparse matrix least square conjugate gradient 

Solution of A'A x = A'B without generate A'A 

inputs:
    A - sparse matrix 
    B - second member of linear system  
    maxIt - max iterations 
    cgTol - tolerance 
*/
float * sparse_lscg(sparseMatrix A, float * B, int maxIt, float cgTol)
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
            q[row] = 0.0f;                           // q = 0, for multiplication

        for (int ind = 0; ind < A.nnz; ind++) 
            q[A.i[ind]] += A.v[ind] * p[A.j[ind]];   // q = G * p   
    }
    
    delete[] s; delete[] q; delete[] r; delete[] p;
    
    return x;
}

/* */
float * movingAverageSmoothing(float * volume, int nx, int ny, int nz, int samples)
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

/* */
float * gaussianFilterSmoothing(float * volume, int nx, int ny, int nz, float stdv, int samples)
{
    int init = samples / 2;
    int nPoints = nx*ny*nz;
    int nKernel = samples*samples*samples;

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

/* */
sparseMatrix getDerivativeMatrix(int n, int degree)
{
	sparseMatrix L;
    
    int elements = degree + 1;    // number of elements per line

    L.init(n, n - degree, elements * n);

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

/* Function to compute a trilinear interpolation */
float triLinearInterpolation(float c000, float c001, float c100, float c101, float c010, float c011, float c110, float c111, 
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

# endif
