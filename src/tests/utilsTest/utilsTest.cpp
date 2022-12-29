# include <string>
# include <vector>
# include <iostream>

# include "../../essentials/utils.hpp"

int main(int argc, char**argv)
{
    // Testing min - max algorithms ----------------------------------------------------------------

    std::cout<<"Testing min - max algorithms\n"<<std::endl;

    std::cout<<"\nMin integer between 5 and 3 is: "<<imin(5, 3)<<std::endl;
    std::cout<<"Min integer between 3 and 5 is: "<<imin(3, 5)<<std::endl;
    
    std::cout<<"\nMax integer between 2 and 6 is: "<<imax(2, 6)<<std::endl;
    std::cout<<"Max integer between 6 and 2 is: "<<imax(6, 2)<<std::endl;

    std::cout<<"\nMin float between 5.3 and 3.5 is: "<<min(5.3f, 3.5f)<<std::endl;
    std::cout<<"Min float between 3.5 and 5.3 is: "<<min(3.5f, 5.3f)<<std::endl;

    std::cout<<"\nMax float between 2.6 and 6.2 is: "<<max(2.6f, 6.2f)<<std::endl;
    std::cout<<"Max float between 6.2 and 2.6 is: "<<max(6.2f, 2.6f)<<std::endl;

    std::cout<<"\nMin float between 5.3, 3.3 and 3.5 is: "<<min3(5.3f, 3.3f, 3.5f)<<std::endl;
    std::cout<<"Max float between 3.5, 3.3 and 5.3 is: "<<max3(3.5f, 3.3f, 5.3f)<<std::endl;

    std::cout<<"\nMin float between 5.5, 5.3, 3.3 and 3.5 is: "<<min4(5.5f, 5.3f, 3.3f, 3.5f)<<std::endl;

    // Testing pick up parameter from file -------------------------------------------------------

    std::cout<<"\nTesting pick up parameter from file"<<std::endl;

    std::string file = "outputs/parametersTest.txt";

    // string reading structure
    std::string p1 = catchParameter("parameter1", file);         

    // integer reading structure
    int p2 = std::stoi(catchParameter("parameter2", file));

    // float reading structure
    float p3 = std::stof(catchParameter("parameter3", file));

    // vector reading structure
    std::vector<std::string> p4 = split(catchParameter("parameter4", file), ',');
    
    std::vector<float> p4f = {};
    for (auto i : p4) 
        p4f.push_back(std::stof(i));

    // boolean test
    bool p5 = str2bool(catchParameter("parameter5", file)); 

    // Printing all tests

    std::cout<<"\nparameter1 = "<<p1<<std::endl;
    std::cout<<"parameter2 = "<<p2<<std::endl;
    std::cout<<"parameter3 = "<<p3<<std::endl;
    std::cout<<"parameter4 = "<<p4f[0]<<", "<<p4f[1]<<", "<<p4f[2]<<std::endl;
    std::cout<<"parameter5 = "<<p5<<"\n"<<std::endl;

    // Testing sparse matrix ------------------------------------------------------------
    std::cout<<"Testing sparse matrix\n"<<std::endl;

    std::cout<<"Objective: Solve the linear system below"<<std::endl;
    std::cout<<"A =         x =     b ="<<std::endl;
    std::cout<<"    2 1 3       1        6"<<std::endl;
    std::cout<<"    5 4 2       1       11"<<std::endl;
    std::cout<<"    4 2 1       1        7"<<std::endl;
    std::cout<<"Target functions:"<<std::endl;
    std::cout<<"    - Struct sparseMatrix"<<std::endl;
    std::cout<<"    - sparse_lscg()"<<std::endl;
    std::cout<<"\nSolution: [1, 1, 1]\n"<<std::endl;

    sparseMatrix A;

    int nrows = 3;
    int ncols = 3;
    int nonZeros = 9;

    A.init(nrows, ncols, nonZeros);

    std::vector<float> values = {2.0f, 1.0f, 3.0f, 5.0f, 4.0f, 2.0f, 4.0f, 2.0f, 1.0f};

    int index = 0;

    for (int i = 0; i < A.n; i++)
    {    
        for(int j = 0; j < A.m; j++)
        {
            A.i[index] = i;
            A.j[index] = j;        
            A.v[index] = values[index];

            std::cout<<"A row = "<<i<<", A col = "<<j<<", A value = "<<values[index]<<"\n";

            index += 1;    
        }
    }

    float * B = new float[3];

    B[0] = 6.0f; B[1] = 11.0f; B[2] = 7.0f;

    float * xp = sparse_lscg(A, B, 10, 1e-10);

    std::cout<<"\nx = "<<xp[0]<<", "<<xp[1]<<", "<<xp[2]<<std::endl;

    // Testing derivative matrix generator ---------------------------------------------------------------------

    std::cout<<"\nTesting derivative matrix generator"<<std::endl;

    std::cout<<"\nL 0th order = identity"<<std::endl;

    sparseMatrix L = getDerivativeMatrix(5, 0);

    for (int index = 0; index < L.nnz; index++)
    {
        std::cout<<"row = "<<L.i[index]<<", col = "<<L.j[index]<<", value = "<<L.v[index]<<"\n";
    }

    std::cout<<"\nL 1st order"<<std::endl;

    L = getDerivativeMatrix(5, 1);
    
    for (int index = 0; index < L.nnz; index++)
    {
        std::cout<<"row = "<<L.i[index]<<", col = "<<L.j[index]<<", value = "<<L.v[index]<<"\n";
    }
    std::cout<<"\nL 2nd order"<<std::endl;

    L = getDerivativeMatrix(5, 2);
    
    for (int index = 0; index < L.nnz; index++)
    {
        std::cout<<"row = "<<L.i[index]<<", col = "<<L.j[index]<<", value = "<<L.v[index]<<"\n";
    }

    // Testing smoothing filters --------------------------------------------------------------------------

    std::cout<<"\nTesting smoothing filters"<<std::endl;

    int nx = 31;
    int ny = 31;
    int nz = 31;

    float * volume = new float[nx*ny*nz]();

    volume[15 + 15*nz + 15*nz*nx] = 100;

    float * mvSmooth = movingAverageSmoothing(volume, nx, ny, nz, 5); 
    float * gnSmooth = gaussianFilterSmoothing(volume, nx, ny, nz, 1.5, 5);

    writeBinaryFloat("outputs/mvSmooth.bin", mvSmooth, nx*ny*nz);
    writeBinaryFloat("outputs/gnSmooth.bin", gnSmooth, nx*ny*nz);

    std::cout<<"\nSmoothed volumes written"<<std::endl;

    // Testing trilinear interpolation ---------------------------------------------------------------------

    std::cout<<"\nTesting trilinear interpolation"<<std::endl;

    std::cout<<"\nF(z,x,y) = 1500 + 90*z + 30*x - 30*y"<<std::endl;

    float z = 0.5f;
    float x = 0.7f;
    float y = 0.2f;

    std::cout<<"\nPoint: (z = 0.5, 0.7, 0.2)"<<std::endl;
    std::cout<<"\nknowing: F(0,0,0), F(0,0,1), F(0,1,0), F(1,0,0)"<<std::endl;
    std::cout<<"          F(1,1,0), F(1,0,1), F(1,1,0), F(1,1,1)"<<std::endl;

    float fo = 1500.0f;
    float fz = 90.0f;
    float fx = 30.0f;
    float fy =-30.0f;

    float f_exact = fo + fz*z + fx*x + fy*y; 

    std::cout<<"\nExact solution = "<<f_exact<<std::endl;

    float x0 = 0.0f; float x1 = 1.0f; 
    float y0 = 0.0f; float y1 = 1.0f; 
    float z0 = 0.0f; float z1 = 1.0f; 

    float c000 = fo + fz*z0 + fx*x0 + fy*y0;
    float c001 = fo + fz*z1 + fx*x0 + fy*y0;
    float c010 = fo + fz*z0 + fx*x0 + fy*y1;
    float c100 = fo + fz*z0 + fx*x1 + fy*y0;
    float c101 = fo + fz*z1 + fx*x1 + fy*y0;
    float c011 = fo + fz*z1 + fx*x0 + fy*y1;
    float c110 = fo + fz*z0 + fx*x0 + fy*y0;
    float c111 = fo + fz*z1 + fx*x1 + fy*y1;

    float f_interp = triLinearInterpolation(c000, c001, c100, c101, c010, c011, c110, c111, x0, x1, y0, y1, z0, z1, x, y, z);

    std::cout<<"Interpolated solution = "<<f_interp<<std::endl;

    return 0;
}
