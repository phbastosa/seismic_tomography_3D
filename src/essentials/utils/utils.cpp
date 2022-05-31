# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <algorithm>

# include "utils.hpp"
# include "../model/model.hpp"

int Utils::imin(int v1, int v2)
{
    if (v1 < v2)
        return v1;
    else
        return v2;
}

int Utils::imax(int v1, int v2)
{
    if (v1 > v2)
        return v1;
    else
        return v2;
}

float Utils::min(float v1, float v2)
{
    if (v1 < v2)
        return v1;
    else
        return v2;
}

float Utils::max(float v1, float v2)
{
    if (v1 > v2)
        return v1;
    else
        return v2;
}

float Utils::max3(float v1, float v2, float v3)
{
    float max = v1;

    if (max < v2) max = v2;
    
    if (max < v3) max = v3;

    return max;
}

float Utils::min3(float v1, float v2, float v3)
{
    float min = v1;

    if (min > v2) min = v2;
    
    if (min > v3) min = v3;

    return min;
}

float Utils::min4(float v1, float v2, float v3, float v4)
{
    float min = v1;
    
    if (min > v2) min = v2;
    if (min > v3) min = v3;
    if (min > v4) min = v4;

    return min;
}

std::vector<float> Utils::linspace(float start, float end, int num)
{
    std::vector<float> linspaced;
    
    if (num == 0) return linspaced;
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    } 

    linspaced.reserve(num);

    float delta = (end - start) / (num - 1);

    for (int i = 0; i < num; i++)
    {
        linspaced.emplace_back(start + (float)(delta*i));
    }

    return linspaced;
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

float * Utils::sparse_cgls(int * iG, int * jG, float * vG, float * B, int nD, int nM, int nnz, int maxIt, float cgTol)
{
    float a, b, qTq, rTr, rd;

    float * s = new float[nD]();
    float * q = new float[nD]();
    float * r = new float[nM]();
    float * p = new float[nM]();
    float * x = new float[nM]();    // Linear system solution

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int i = 0; i < nD; i++) 
        s[i] = B[i]; 

    // r = G' * s    
    for (int i = 0; i < nnz; i++) 
        r[jG[i]] += vG[i] * s[iG[i]];        

    // p = r
    for (int i = 0; i < nM; i++) 
        p[i] = r[i]; 

    // q = G * p
    for (int i = 0; i < nnz; i++) 
        q[iG[i]] += vG[i] * p[jG[i]];        

    for (int i = 0; i < maxIt; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < nD; k++)          // q inner product
            qTq += q[k] * q[k];               // qTq = q' * q

        rTr = 0.0f;
        for (int k = 0; k < nM; k++)          // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r 

        a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

        for (int k = 0; k < nM; k++)          // model atualization
            x[k] += a * p[k];                 // x = x + a * p

        for (int k = 0; k < nD; k++)          // s atualization  
            s[k] -= a * q[k];                 // s = s - a * q 

        rd = 0.0f;
        for (int k = 0; k < nM; k++)          // r inner product for division 
            rd += r[k] * r[k];                // rd = r' * r

        for (int k = 0; k < nM; k++)          // Zeroing r 
            r[k] = 0.0f;                      // r = 0, for multiplication
        
        for (int k = 0; k < nnz; k++)         // r atualization 
            r[jG[k]] += vG[k] * s[iG[k]];     // r = G' * s    

        rTr = 0.0f;                
        for (int k = 0; k < nM; k++)          // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r

        if (sqrtf(rd) < cgTol) break;          // Convergence condition
        
        b = rTr / rd;                         // b = (r' * r) / rd

        for (int k = 0; k < nM; k++)          
            p[k] = r[k] + b * p[k];           // p = r + b * p 

        for (int k = 0; k < nD; k++) 
            q[k] = 0.0f;                      // q = 0, for multiplication

        for (int k = 0; k < nnz; k++) 
            q[iG[k]] += vG[k] * p[jG[k]];     // q = G * p   
    }
    
    delete[] s; delete[] q; delete[] r; delete[] p;
    
    return x;
}

float Utils::triLinearInterpolation(Utils::point3D p3D, Model m3D, float *T)
{
    float x0 = floorf(p3D.x / m3D.dx) * m3D.dx;
    float y0 = floorf(p3D.y / m3D.dy) * m3D.dy;
    float z0 = floorf(p3D.z / m3D.dz) * m3D.dz;

    float x1 = floorf(p3D.x / m3D.dx) * m3D.dx + m3D.dx;
    float y1 = floorf(p3D.y / m3D.dy) * m3D.dy + m3D.dy;
    float z1 = floorf(p3D.z / m3D.dz) * m3D.dz + m3D.dz;

    int xi = (int)(p3D.x / m3D.dx) + m3D.nb;    
    int yi = (int)(p3D.y / m3D.dy) + m3D.nb;    
    int zi = (int)(p3D.z / m3D.dz) + m3D.nb;    

    int indT = zi + xi*m3D.nzz + yi*m3D.nxx*m3D.nzz;

    float c000 = 1.0f / T[indT];
    float c001 = 1.0f / T[indT + 1];
    float c100 = 1.0f / T[indT + m3D.nzz]; 
    float c101 = 1.0f / T[indT + 1 + m3D.nzz]; 
    float c010 = 1.0f / T[indT + m3D.nxx*m3D.nzz]; 
    float c011 = 1.0f / T[indT + 1 + m3D.nxx*m3D.nzz]; 
    float c110 = 1.0f / T[indT + m3D.nzz + m3D.nxx*m3D.nzz]; 
    float c111 = 1.0f / T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz];

    float xd = (p3D.x - x0) / (x1 - x0);
    float yd = (p3D.y - y0) / (y1 - y0);
    float zd = (p3D.z - z0) / (z1 - z0);

    float c00 = c000*(1 - xd) + c100*xd;    
    float c01 = c001*(1 - xd) + c101*xd;    
    float c10 = c010*(1 - xd) + c110*xd;    
    float c11 = c011*(1 - xd) + c111*xd;    

    float c0 = c00*(1 - yd) + c10*yd;
    float c1 = c01*(1 - yd) + c11*yd;

    return 1.0f / (c0*(1 - zd) + c1*zd);
}

float Utils::triCubicInterpolation(Utils::point3D p3D, Model m3D, float *T)
{
    std::vector< int > iM, jM;
    std::vector<float> vM, b;

    int n = 64;
    int nnz = 1000;
    float result = 0;

    b.reserve(n);

    iM = {0,1,2,2,2,2,3,3,3,3,4,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
         16,17,18,18,18,18,19,19,19,19,20,21,22,22,22,22,23,23,23,23,24,24,24,24,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,28,28,28,28,29,29,29,29,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
         32,32,32,32,33,33,33,33,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,34,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,36,36,36,36,37,37,37,37,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,41,41,41,41,
         41,41,41,41,41,41,41,41,41,41,41,41,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,42,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,
         43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,
         46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,
         48,48,48,48,49,49,49,49,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51,52,52,52,52,53,53,53,53,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,56,57,57,57,57,
         57,57,57,57,57,57,57,57,57,57,57,57,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,58,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,
         59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,61,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,
         62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63};

    jM = {0,8,0,1,8,9,0,1,8,9,16,32,16,17,32,33,16,17,32,33,0,2,16,18,8,10,32,34,0,1,2,3,8,9,10,11,16,17,18,19,32,33,34,35,0,1,2,3,8,9,10,11,16,17,18,19,32,33,34,35,0,2,16,18,8,10,32,34,0,1,2,3,8,9,10,11,16,17,18,19,32,33,34,35,0,1,2,3,8,9,10,11,16,17,18,19,32,33,34,35,
         24,40,24,25,40,41,24,25,40,41,48,56,48,49,56,57,48,49,56,57,24,26,48,50,40,42,56,58,24,25,26,27,40,41,42,43,48,49,50,51,56,57,58,59,24,25,26,27,40,41,42,43,48,49,50,51,56,57,58,59,24,26,48,50,40,42,56,58,24,25,26,27,40,41,42,43,48,49,50,51,56,57,58,59,24,25,26,27,40,41,42,43,48,49,50,51,56,57,58,59,
          0,4,24,28,8,12,40,44,0,1,4,5,8,9,12,13,24,25,28,29,40,41,44,45,0,1,4,5,8,9,12,13,24,25,28,29,40,41,44,45,16,20,48,52,32,36,56,60,16,17,20,21,32,33,36,37,48,49,52,53,56,57,60,61,16,17,20,21,32,33,36,37,48,49,52,53,56,57,60,61,0,2,4,6,16,18,20,22,24,26,28,30,48,50,52,54,8,10,12,14,
         32,34,36,38,40,42,44,46,56,58,60,62,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
         24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,2,4,6,16,18,20,22,24,26,28,30,48,50,52,54,8,10,12,14,32,34,36,38,40,42,44,46,56,58,60,62,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
         28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
          0,4,24,28,8,12,40,44,0,1,4,5,8,9,12,13,24,25,28,29,40,41,44,45,0,1,4,5,8,9,12,13,24,25,28,29,40,41,44,45,16,20,48,52,32,36,56,60,16,17,20,21,32,33,36,37,48,49,52,53,56,57,60,61,16,17,20,21,32,33,36,37,48,49,52,53,56,57,60,61,0,2,4,6,16,18,20,22,24,26,28,30,48,50,52,54,8,10,12,14,
         32,34,36,38,40,42,44,46,56,58,60,62,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
         24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,2,4,6,16,18,20,22,24,26,28,30,48,50,52,54,8,10,12,14,32,34,36,38,40,42,44,46,56,58,60,62,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
         28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};

    vM = {1,1,-3,3,-2,-1,2,-2,1,1,1,1,-3,3,-2,-1,2,-2,1,1,-3,3,-2,-1,-3,3,-2,-1,9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,2,-2,1,1,2,-2,1,1,-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,
          1,1,-3,3,-2,-1,2,-2,1,1,1,1,-3,3,-2,-1,2,-2,1,1,-3,3,-2,-1,-3,3,-2,-1,9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,2,-2,1,1,2,-2,1,1,-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,
         -3,3,-2,-1,-3,3,-2,-1,9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,-3,3,-2,-1,-3,3,-2,-1,9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,9,-9,-9,9,
          6,3,-6,-3,6,-6,3,-3,4,2,2,1,-27,27,27,-27,27,-27,-27,27,-18,-9,18,9,18,9,-18,-9,-18,18,-9,9,18,-18,9,-9,-18,18,18,-18,-9,9,9,-9,-12,-6,-6,-3,12,6,6,3,-12,-6,12,6,-6,-3,6,3,-12,12,-6,6,-6,6,-3,3,-8,-4,-4,-2,-4,-2,-2,-1,18,-18,-18,18,-18,18,18,-18,9,9,-9,-9,-9,-9,9,9,12,-12,6,-6,-12,12,-6,6,
         12,-12,-12,12,6,-6,-6,6,6,6,3,3,-6,-6,-3,-3,6,6,-6,-6,3,3,-3,-3,8,-8,4,-4,4,-4,2,-2,4,4,2,2,2,2,1,1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,18,-18,-18,18,-18,18,18,-18,12,6,-12,-6,-12,-6,12,6,9,-9,9,-9,-9,9,-9,9,12,-12,-12,12,
          6,-6,-6,6,6,3,6,3,-6,-3,-6,-3,8,4,-8,-4,4,2,-4,-2,6,-6,6,-6,3,-3,3,-3,4,2,4,2,2,1,2,1,-12,12,12,-12,12,-12,-12,12,-6,-6,6,6,6,6,-6,-6,-6,6,-6,6,6,-6,6,-6,-8,8,8,-8,-4,4,4,-4,-3,-3,-3,-3,3,3,3,3,-4,-4,4,4,-2,-2,2,2,-4,4,-4,4,-2,2,-2,2,-2,-2,-2,-2,-1,-1,-1,-1,
          2,-2,1,1,2,-2,1,1,-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,2,-2,1,1,2,-2,1,1,-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,-6,6,6,-6,
         -4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,18,-18,-18,18,-18,18,18,-18,12,6,-12,-6,-12,-6,12,6,12,-12,6,-6,-12,12,-6,6,9,-9,-9,9,9,-9,-9,9,8,4,4,2,-8,-4,-4,-2,6,3,-6,-3,6,3,-6,-3,6,-6,3,-3,6,-6,3,-3,4,2,2,1,4,2,2,1,-12,12,12,-12,12,-12,-12,12,-6,-6,6,6,6,6,-6,-6,-8,8,-4,4,8,-8,4,-4,
         -6,6,6,-6,-6,6,6,-6,-4,-4,-2,-2,4,4,2,2,-3,-3,3,3,-3,-3,3,3,-4,4,-2,2,-4,4,-2,2,-2,-2,-1,-1,-2,-2,-1,-1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1,-12,12,12,-12,12,-12,-12,12,-8,-4,8,4,8,4,-8,-4,-6,6,-6,6,6,-6,6,-6,-6,6,6,-6,
         -6,6,6,-6,-4,-2,-4,-2,4,2,4,2,-4,-2,4,2,-4,-2,4,2,-3,3,-3,3,-3,3,-3,3,-2,-1,-2,-1,-2,-1,-2,-1,8,-8,-8,8,-8,8,8,-8,4,4,-4,-4,-4,-4,4,4,4,-4,4,-4,-4,4,-4,4,4,-4,-4,4,4,-4,-4,4,2,2,2,2,-2,-2,-2,-2,2,2,-2,-2,2,2,-2,-2,2,-2,2,-2,2,-2,2,-2,1,1,1,1,1,1,1,1};

    int xi = (int)(p3D.x / m3D.dx) + m3D.nb;    
    int yi = (int)(p3D.y / m3D.dy) + m3D.nb;    
    int zi = (int)(p3D.z / m3D.dz) + m3D.nb;    

    int indT = zi + xi*m3D.nzz + yi*m3D.nxx*m3D.nzz;

    // vs of f(x,y,z) at each corner
    b.emplace_back(T[indT]);
    b.emplace_back(T[indT + m3D.nzz]);
    b.emplace_back(T[indT + m3D.nxx*m3D.nzz]);
    b.emplace_back(T[indT + m3D.nzz + m3D.nxx*m3D.nzz]);
    b.emplace_back(T[indT + 1]);
    b.emplace_back(T[indT + 1 + m3D.nzz]);
    b.emplace_back(T[indT + 1 + m3D.nxx*m3D.nzz]);
    b.emplace_back(T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz]);
    // vs of df/dx at each corner
    b.emplace_back(0.5f * (T[indT + m3D.nzz] - T[indT - m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2*m3D.nzz] - T[indT]));
    b.emplace_back(0.5f * (T[indT + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT - m3D.nzz + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz] - T[indT + 1 - m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + 2*m3D.nzz] - T[indT + 1]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nxx*m3D.nzz]));
    // vs of df/dy at each corner
    b.emplace_back(0.5f * (T[indT + m3D.nxx*m3D.nzz] - T[indT - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2*m3D.nxx*m3D.nzz] - T[indT]));
    b.emplace_back(0.5f * (T[indT + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + 2*m3D.nxx*m3D.nzz] - T[indT + 1]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz]));
    // vs of df/dz at each corner.
    b.emplace_back(0.5f * (T[indT + 1] - T[indT - 1]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz] - T[indT - 1 + m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nzz + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2] - T[indT]));
    b.emplace_back(0.5f * (T[indT + 2 + m3D.nzz] - T[indT + m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2 + m3D.nxx*m3D.nzz] - T[indT + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.5f * (T[indT + 2 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nzz + m3D.nxx*m3D.nzz]));
    // vs of d2f/dxdy at each corner
    b.emplace_back(0.25f * (T[indT + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nzz - m3D.nxx*m3D.nzz] + T[indT - m3D.nzz - m3D.nxx*m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nxx*m3D.nzz] - T[indT + 2*m3D.nzz - m3D.nxx*m3D.nzz] + T[indT - m3D.nxx*m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT - m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + m3D.nzz] + T[indT - m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 2*m3D.nxx*m3D.nzz] - T[indT + 2*m3D.nzz] + T[indT]));    
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 1 - m3D.nzz - m3D.nxx*m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 1 - m3D.nxx*m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz] + T[indT + 1 - m3D.nzz]));    
    b.emplace_back(0.25f * (T[indT + 1 + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nzz] + T[indT + 1]));    
    // vs of d2f/dxdz at each corner
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz] - T[indT + 1 - m3D.nzz] - T[indT - 1 + m3D.nzz] + T[indT - 1 - m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 1 + 2*m3D.nzz] - T[indT + 1] - T[indT - 1 + 2*m3D.nzz] + T[indT - 1]));
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 - m3D.nzz + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT - 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + m3D.nzz] - T[indT + 2 - m3D.nzz] - T[indT + m3D.nzz] + T[indT - m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + 2*m3D.nzz] - T[indT + 2] - T[indT + 2*m3D.nzz] + T[indT]));
    b.emplace_back(0.25f * (T[indT + 2 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - m3D.nzz + m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nxx*m3D.nzz] - T[indT + 2*m3D.nzz + m3D.nxx*m3D.nzz] + T[indT + m3D.nxx*m3D.nzz]));
    // vs of d2f/dydz at each corner
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nxx*m3D.nzz] + T[indT - 1 - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 1 + 2*m3D.nxx*m3D.nzz] - T[indT + 1] - T[indT - 1 + 2*m3D.nxx*m3D.nzz] + T[indT - 1]));
    b.emplace_back(0.25f * (T[indT + 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz] - T[indT - 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + m3D.nxx*m3D.nzz] - T[indT + 2 - m3D.nxx*m3D.nzz] - T[indT + m3D.nxx*m3D.nzz] + T[indT - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nzz - m3D.nxx*m3D.nzz] - T[indT + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT + m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.25f * (T[indT + 2 + 2*m3D.nxx*m3D.nzz] - T[indT + 2] - T[indT + 2*m3D.nxx*m3D.nzz] + T[indT]));
    b.emplace_back(0.25f * (T[indT + 2 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nzz] - T[indT + m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT + m3D.nzz]));
    // vs of d3f/dxdydz at each corner
    b.emplace_back(0.125f * (T[indT + 1 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 1 - m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - 1 + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 - m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - 1 - m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 1 - m3D.nxx*m3D.nzz] - T[indT - 1 + 2*m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nxx*m3D.nzz] + T[indT - 1 + 2*m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - 1 - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 - m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + m3D.nzz] + T[indT + 1 - m3D.nzz] - T[indT - 1 + m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT - 1 - m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT - 1 + m3D.nzz] - T[indT - 1 - m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 1 + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nxx*m3D.nzz] - T[indT + 1 + 2*m3D.nzz] + T[indT + 1] - T[indT - 1 + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT - 1 + 2*m3D.nxx*m3D.nzz] + T[indT - 1 + 2*m3D.nzz] - T[indT - 1]));
    b.emplace_back(0.125f * (T[indT + 2 + m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 - m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 2 - m3D.nzz - m3D.nxx*m3D.nzz] - T[indT + m3D.nzz + m3D.nxx*m3D.nzz] + T[indT - m3D.nzz + m3D.nxx*m3D.nzz] + T[indT + m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - m3D.nzz - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 2 + 2*m3D.nzz + m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nxx*m3D.nzz] - T[indT + 2 + 2*m3D.nzz - m3D.nxx*m3D.nzz] + T[indT + 2 - m3D.nxx*m3D.nzz] - T[indT + 2*m3D.nzz + m3D.nxx*m3D.nzz] + T[indT + m3D.nxx*m3D.nzz] + T[indT + 2*m3D.nzz - m3D.nxx*m3D.nzz] - T[indT - m3D.nxx*m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 2 + m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 2 - m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 2 + m3D.nzz] + T[indT + 2 - m3D.nzz] - T[indT + m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT - m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT + m3D.nzz] - T[indT - m3D.nzz]));
    b.emplace_back(0.125f * (T[indT + 2 + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] - T[indT + 2 + 2*m3D.nxx*m3D.nzz] - T[indT + 2 + 2*m3D.nzz] + T[indT + 2] - T[indT + 2*m3D.nzz + 2*m3D.nxx*m3D.nzz] + T[indT + 2*m3D.nxx*m3D.nzz] + T[indT + 2*m3D.nzz] - T[indT]));

    float * alpha = new float[n](); 

    for (int i = 0; i < nnz; i++) alpha[iM[i]] += vM[i] * b[jM[i]];

    for (int k = 0; k < 4; k++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int i = 0; i < 4; i++)
            {
                int ijk = i + 4*j + 16*k;

                result += alpha[ijk] * powf(p3D.x, i) * powf(p3D.y, j) * powf(p3D.z, k);    
            }
        } 
    }

    delete[] alpha;

    return result;
}
