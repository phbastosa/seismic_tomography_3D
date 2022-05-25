# include <cmath>
# include <vector>
# include <string>
# include <sstream>
# include <algorithm>

# include "utils.hpp"

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

