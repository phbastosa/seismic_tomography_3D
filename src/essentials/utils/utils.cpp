# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <algorithm>

# include "utils.hpp"
# include "../model/model.hpp"

int Utils::imin(int v1, int v2) { return !(v1 > v2) ? v1 : v2; }

int Utils::imax(int v1, int v2) { return !(v1 < v2) ? v1 : v2; }

float Utils::min(float v1, float v2) { return !(v1 > v2) ? v1 : v2; }

float Utils::max(float v1, float v2) { return !(v1 < v2) ? v1 : v2; }

float Utils::max3(float v1, float v2, float v3) { return max(v1, max(v2, v3)); }

float Utils::min3(float v1, float v2, float v3) { return min(v1, min(v2, v3)); }

float Utils::min4(float v1, float v2, float v3, float v4) { return min(v1, min(v2, min(v3, v4))); }

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

float Utils::triLinearInterpolation(Point p, Model m, float * volume)
{
    float x0 = floorf(p.x / m.dx) * m.dx;
    float y0 = floorf(p.y / m.dy) * m.dy;
    float z0 = floorf(p.z / m.dz) * m.dz;

    float x1 = floorf(p.x / m.dx) * m.dx + m.dx;
    float y1 = floorf(p.y / m.dy) * m.dy + m.dy;
    float z1 = floorf(p.z / m.dz) * m.dz + m.dz;

    int xi = (int)(p.x / m.dx) + m.nb;    
    int yi = (int)(p.y / m.dy) + m.nb;    
    int zi = (int)(p.z / m.dz) + m.nb;    

    int id = zi + xi*m.nzz + yi*m.nxx*m.nzz;

    float c000 = 1.0f / volume[id];
    float c001 = 1.0f / volume[id + 1];
    float c100 = 1.0f / volume[id + m.nzz]; 
    float c101 = 1.0f / volume[id + 1 + m.nzz]; 
    float c010 = 1.0f / volume[id + m.nxx*m.nzz]; 
    float c011 = 1.0f / volume[id + 1 + m.nxx*m.nzz]; 
    float c110 = 1.0f / volume[id + m.nzz + m.nxx*m.nzz]; 
    float c111 = 1.0f / volume[id + 1 + m.nzz + m.nxx*m.nzz];

    float xd = (p.x - x0) / (x1 - x0);
    float yd = (p.y - y0) / (y1 - y0);
    float zd = (p.z - z0) / (z1 - z0);

    float c00 = c000*(1 - xd) + c100*xd;    
    float c01 = c001*(1 - xd) + c101*xd;    
    float c10 = c010*(1 - xd) + c110*xd;    
    float c11 = c011*(1 - xd) + c111*xd;    

    float c0 = c00*(1 - yd) + c10*yd;
    float c1 = c01*(1 - yd) + c11*yd;

    return 1.0f / (c0*(1 - zd) + c1*zd);
}
