# include <cmath>

# include "utils.hpp"

# include "../inout/inout.hpp"

int main(int argc, char **argv)
{
    auto utils = Utils();
    
    std::vector<float> x = utils.linspace(0.25f,0.75f,5);

    for (int i = 0; i < x.size(); i++) std::cout<<x[i]<<" ";

    Utils::Point p;

    p.x = 10.0f;
    p.y = 12.4f;
    p.z = 19.8f;

    std::cout<<"\n\nPoint 3D ("<<p.x<<", "<<p.y<<" ,"<<p.z<<") correctly written!"<<std::endl;

    std::cout<<"\nmax of 2 and 5 is "<<utils.imax(2,5)<<std::endl;
    std::cout<<"min of 2 and 5 is "<<utils.imin(2,5)<<std::endl;
    std::cout<<"max of 2.0f and 5.0f is "<<utils.max(2.0f,5.0f)<<std::endl;
    std::cout<<"min of 2.0f and 5.0f is "<<utils.min(2.0f,5.0f)<<std::endl;

    // Little linear tomography

    float h = 1000.0f;

    float d1 = h;
    float d2 = h * sqrtf(2.0f);
    float d3 = sqrtf((0.5f*h)*(0.5f*h) + h*h);

    int nnz = 57;
    int nD = 19;
    int nM = 9;

    std::vector<int> iA {0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,
                       10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,
                       16,17,17,17,18,18,18};

    std::vector<int> jA {0,1,5,0,4,8,0,1,2,0,1,2,0,4,5,0,1,2,3,4,8,3,4,5,1,2,3,3,4,
                        5,3,7,8,3,4,5,2,3,4,6,7,8,4,5,6,6,7,8,6,7,8,2,4,6,5,6,7};    

    std::vector<float> vA {d3,d3,d3,d2,d2,d2,d1,d1,d1,d1,d1,d1,d3,d3,d3,d1,d1,d1,d3,
                          d3,d3,d1,d1,d1,d3,d3,d3,d1,d1,d1,d3,d3,d3,d1,d1,d1,d3,d3,
                          d3,d1,d1,d1,d3,d3,d3,d1,d1,d1,d1,d1,d1,d2,d2,d2,d3,d3,d3};

    std::vector<float> t {2.189f, 2.612f, 2.000f, 2.000f, 2.142f, 2.000f, 2.018f, 
                          1.875f, 2.189f, 1.875f, 1.941f, 1.875f, 2.142f, 1.666f, 
                          2.018f, 1.666f, 1.666f, 2.612f, 1.941f};

    int * iG = new int[nnz]();
    int * jG = new int[nnz]();
    float * vG = new float[nnz]();
    float * B = new float[nD];

    for (int i = 0; i < nnz; i++)
    {
        iG[i] = iA[i];
        jG[i] = jA[i];
        vG[i] = vA[i];
    }

    for (int i = 0; i < nD; i++) B[i] = t[i];

    int maxIt = 10;
    int cgTol = 1e-6;

    std::vector<float> xTrue {1500, 1500, 1500, 1600, 1600, 1600, 1800, 1800, 1800};

    float * xCalc = Utils::sparse_cgls(iG,jG,vG,B,nD,nM,nnz,maxIt,cgTol);

    std::cout<<"\n\nTrue model!"<<std::endl;
    for (int i = 0; i < nM; i++)    
        std::cout<<xTrue[i]<<" ";
    std::cout<<"\n";

    std::cout<<"\nEstimated model!"<<std::endl;
    for (int i = 0; i < nM; i++)    
        std::cout<<1.0f/xCalc[i]<<" ";
    std::cout<<"\n";

    std::cout<<"\nModel difference!"<<std::endl;
    for (int i = 0; i < nM; i++)    
        std::cout<<(1.0f/xCalc[i]) - xTrue[i]<<" ";
    std::cout<<"\n";

    return 0;
}
