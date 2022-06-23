# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>
# include <string>

class Utils
{
public:
    typedef struct 
    {
        int * i;       // Rows indexes
        int * j;       // Cols indexes 
        float * v;     // Value

        int n;         // Rows number
        int m;         // Cols number
        int nnz;       // Non-zero elements
    
    } sparseMatrix;

    /* */
    float * readBinaryFloat(std::string path, int n);

    /* */
    float * rickerGeneration(int ns, float dt, float fmax);

    /* */
    float * intRickerGeneration(int ns, float dt, float fmax);

    /* */
    void writeBinaryFloat(std::string path, float *array, int n);

    /* */
    std::string catchParameter(std::string target, std::string file);

    /* Function to calculate minimum value between two float inputs */
    float min(float v1, float v2);

    /* Function to calculate maximum value between two float inputs */
    float max(float v1, float v2);

    /* Function to calculate minimum value between two integer inputs */
    int imin(int v1, int v2);
    
    /* Function to calculate maximum value between two integer inputs */
    int imax(int v1, int v2);
       
    /* Function to calculate minimum value between three float inputs */
    float min3(float v1, float v2, float v3);

    /* Function to calculate maximum value between three float inputs */
    float max3(float v1, float v2, float v3);

    /* Function to calculate minimum value between four float inputs */
    float min4(float v1, float v2, float v3, float v4);

    /* Function to convert string to boolean */
    bool str2bool(std::string s);

    /* Function to separete values with a delimiter */
    std::vector<std::string> split(std::string s, char delimiter);

    /* Function to compute a sparse matrix least square conjugate gradient 
    
    Solution of A'A x = A'B without generate A'A 

    inputs:
        A - sparse matrix 
        B - second member of linear system  
        maxIt - max iterations 
        cgTol - tolerance 
    */
    float * sparse_lscg(sparseMatrix A, float * B, int maxIt, float cgTol);

    /* */
    sparseMatrix getDerivativeMatrix(int n, int degree);

    /* Function to compute a trilinear interpolation */
    float triLinearInterpolation(float c000, float c001, float c100, float c101, float c010, float c011, float c110, float c111, 
                                 float x0, float x1, float y0, float y1, float z0, float z1, float x, float y, float z);
};

# endif
