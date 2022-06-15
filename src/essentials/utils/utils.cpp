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

void Utils::readBinaryFloat(std::string path, float *array, int n)
{
    std::ifstream file(path, std::ios::in);
    
    if (file.is_open()) 
    {    
        file.read((char *) array, n * sizeof(float));
    }
    else
    {
        throw std::invalid_argument("Error: file could not be opened!");
    }

    file.close();    
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

float * Utils::sparse_cgls(int * iA, int * jA, float * vA, float * B, int n, int m, int nnz, int maxIt, float cgTol)
{    
    float a, b, qTq, rTr, rd;

    float * s = new float[n]();
    float * q = new float[n]();
    float * r = new float[n]();
    float * p = new float[m]();
    float * x = new float[m]();    // Linear system solution

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int i = 0; i < n; i++) 
        s[i] = B[i]; 

    // r = G' * s    
    for (int i = 0; i < nnz; i++) 
        r[jA[i]] += vA[i] * s[iA[i]];        

    // p = r
    for (int i = 0; i < m; i++) 
        p[i] = r[i]; 

    // q = G * p
    for (int i = 0; i < nnz; i++) 
        q[iA[i]] += vA[i] * p[jA[i]];        

    for (int i = 0; i < maxIt; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < n; k++)          // q inner product
            qTq += q[k] * q[k];               // qTq = q' * q

        rTr = 0.0f;
        for (int k = 0; k < m; k++)          // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r 

        a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

        for (int k = 0; k < m; k++)          // model atualization
            x[k] += a * p[k];                 // x = x + a * p

        for (int k = 0; k < n; k++)          // s atualization  
            s[k] -= a * q[k];                 // s = s - a * q 

        rd = 0.0f;
        for (int k = 0; k < m; k++)          // r inner product for division 
            rd += r[k] * r[k];                // rd = r' * r

        for (int k = 0; k < m; k++)          // Zeroing r 
            r[k] = 0.0f;                      // r = 0, for multiplication
        
        for (int k = 0; k < nnz; k++)         // r atualization 
            r[jA[k]] += vA[k] * s[iA[k]];     // r = G' * s    

        rTr = 0.0f;                
        for (int k = 0; k < m; k++)          // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r

        if (sqrtf(rd) < cgTol) break;          // Convergence condition
        
        b = rTr / rd;                         // b = (r' * r) / rd

        for (int k = 0; k < m; k++)          
            p[k] = r[k] + b * p[k];           // p = r + b * p 

        for (int k = 0; k < n; k++) 
            q[k] = 0.0f;                      // q = 0, for multiplication

        for (int k = 0; k < nnz; k++) 
            q[iA[k]] += vA[k] * p[jA[k]];     // q = G * p   
    }
    
    delete[] s; delete[] q; delete[] r; delete[] p;
    
    return x;
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
