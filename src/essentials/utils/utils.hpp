# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>
# include <string>

# include "../model/model.hpp"

class Utils
{
private:

public:
    typedef struct   // To storage a simple point in space
    {
        float x;     // x position of point  
        float y;     // y position of point
        float z;     // z position of point   
        
    } Point;

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
    static bool str2bool(std::string s);

    /* Function to calculate a linear spaced position */
    static std::vector<float> linspace(float start, float end, int num);        

    /* Function to separete values with a delimiter */
    static std::vector<std::string> split(std::string s, char delimiter);

    /* Function to compute a sparse matrix least square conjugate gradient 
    
    Solution of A'A x = A'B without generate A'A 

    inputs:
        iA - line of non zero elements
        jA - columns of non zero elements
        vA - value of non zero elements
        B - second member of linear system  
        n - number of lines in A
        m - number of columns in A
        nnz - number of non zero elements
        maxIt - max iterations 
        cgTol - tolerance 
    */
    static float * sparse_cgls(int * iA, int * jA, float * vA, float * B, int n, int m, int nnz, int maxIt, float cgTol);

    /* Function to compute a trilinear interpolation */
    static float triLinearInterpolation(Point p, Model m, float * volume);
};

# endif
