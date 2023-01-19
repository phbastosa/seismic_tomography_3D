# include <cmath>
# include <chrono>
# include <vector>
# include <iostream>

# include "../../src/cgls/cgls.cuh"
# include "../../src/essentials/utils.hpp"

int main(int argc, char** argv)
{   
    std::cout<<"Testing matrix operation\n"<<std::endl;

    std::cout<<"Objective: Solve the linear system below"<<std::endl;
    std::cout<<"A =         x =     b ="<<std::endl;
    std::cout<<"    2 1 3       1        6"<<std::endl;
    std::cout<<"    5 4 2       1       11"<<std::endl;
    std::cout<<"    4 2 1       1        7"<<std::endl;
    std::cout<<"\nTarget functions:"<<std::endl;
    std::cout<<"    - Struct sparseMatrix"<<std::endl;
    std::cout<<"    - sparse_cgls_cpu()"<<std::endl;
    std::cout<<"    - sparse_cgls_gpu()"<<std::endl;
    std::cout<<"\nSolution: [1, 1, 1]\n"<<std::endl;

    Utils::sparseMatrix A;

    A.n = 3;
    A.m = 3;
    A.nnz = 9;

    A.i = new int[A.nnz]();
    A.j = new int[A.nnz]();
    A.v = new float[A.nnz]();

    std::vector<float> values = {2.0f, 1.0f, 3.0f, 5.0f, 4.0f, 2.0f, 4.0f, 2.0f, 1.0f};

    int index = 0;

    for (int i = 0; i < 3; i++)
    {    
        for(int j = 0; j < 3; j++)
        {
            A.i[index] = i;
            A.j[index] = j;        
            A.v[index] = values[index];

            std::cout<<"A row = "<<i<<", A col = "<<j<<", A value = "<<values[index]<<"\n";

            index += 1;    
        }
    }

    float * B = new float[A.n];

    B[0] = 6.0f; B[1] = 11.0f; B[2] = 7.0f;

    float * x_cpu = new float[A.m];
    float * x_gpu = new float[A.m];

    sparse_cgls_cpu(A.i, A.j, A.v, B, x_cpu, A.n, A.m, A.nnz, 10, 1e-10);
    sparse_cgls_gpu(A.i, A.j, A.v, B, x_gpu, A.n, A.m, A.nnz, 10, 1e-10);

    std::cout<<"\n Running on CPU: x = "<<x_cpu[0]<<", "<<x_cpu[1]<<", "<<x_cpu[2]<<std::endl;
    std::cout<<" Running on GPU: x = "<<x_cpu[0]<<", "<<x_cpu[1]<<", "<<x_cpu[2]<<std::endl;

    // // Linear tomography example 

    // int M = 9;
    // int N = 19;
    // int NNZ = 57;

    // float h = 1000.0f;

    // float d1 = h;
    // float d2 = h * sqrt(2.0f);
    // float d3 = sqrt((0.5f*h)*(0.5f*h) + h*h);

    // int iM[] = {0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,
    // 	        9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,
    // 	       15,15,16,16,16,17,17,17,18,18,18};

    // int jM[] = {0,1,5,0,4,8,0,1,2,0,1,2,0,4,5,0,1,2,3,4,8,3,4,5,1,2,3,
    //             3,4,5,3,7,8,3,4,5,2,3,4,6,7,8,4,5,6,6,7,8,6,7,8,2,4,6,
    //             5,6,7};    

    // float vM[] = {d3,d3,d3,d2,d2,d2,d1,d1,d1,d1,d1,d1,d3,d3,d3,d1,d1,
    // 			   d1,d3,d3,d3,d1,d1,d1,d3,d3,d3,d1,d1,d1,d3,d3,d3,d1,
    // 			   d1,d1,d3,d3,d3,d1,d1,d1,d3,d3,d3,d1,d1,d1,d1,d1,d1,
    // 			   d2,d2,d2,d3,d3,d3};

    // float b[] = {2.189f,2.612f,2.000f,2.000f,2.142f,2.000f,2.018f, 
    //              1.875f,2.189f,1.875f,1.941f,1.875f,2.142f,1.666f,
    //              2.018f,1.666f,1.666f,2.612f,1.941f};

    // float * x = new float[M]();
    // float * B = new float[N]();
    
    // int * iA = new int[NNZ]();
    // int * jA = new int[NNZ]();
    // float * vA = new float[NNZ]();

    // for (int index = 0; index < NNZ; index++)
    // {
    //     iA[index] = iM[index];
    //     jA[index] = jM[index];
    //     vA[index] = vM[index];
    // }

    // for (int row = 0; row < N; row++)
    //     B[row] = b[row];

    // int NIT = 10;
    // float TOL = 1e-6f;    

    // auto ti = std::chrono::system_clock::now();
    // sparse_cgls_gpu(iA, jA, vA, B, x, N, M, NNZ, NIT, TOL);    
    // auto tf = std::chrono::system_clock::now();

    // std::chrono::duration<double> runTime = tf - ti;

    // std::cout<<"\n"<<std::endl;

    // for (int i = 0; i < M; i++)
    // {
    //     std::cout<<1.0f / x[i]<<std::endl;
    // }

    // std::cout<<"\nRuntime: "<<runTime.count()<<" s."<<std::endl;

    return 0;
}






