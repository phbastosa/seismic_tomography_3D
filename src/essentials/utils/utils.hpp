# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>
# include <string>

# include "../model/model.hpp"

class Utils
{
private:

public:
    
    typedef struct 
    {
        float x;
        float y; 
        
    } point2D;

    typedef struct 
    {
        float x;
        float y; 
        float z; 
        
    } point3D;

    /* Function to calculate minimum value between two inputs */
    float min(float v1, float v2);

    /* Function to calculate maximum value between two inputs */
    float max(float v1, float v2);

    /* */
    int imin(int v1, int v2);
    
    /* */
    int imax(int v1, int v2);
       
    /* */
    float min3(float v1, float v2, float v3);

    /* */
    float max3(float v1, float v2, float v3);

    /* Function to calculate minimum value between four inputs */
    float min4(float v1, float v2, float v3, float v4);

    static bool str2bool(std::string s);

    static std::vector<float> linspace(float start, float end, int num);        

    static std::vector<std::string> split(std::string s, char delimiter);

    static float * sparse_cgls(int * iG, int * jG, float * vG, float * B, int nD, int nM, int nnz, int maxIt, float cgTol);

    static float triCubicInterpolation(Utils::point3D p3D, Model m3D, float *T);

    static float triLinearInterpolation(Utils::point3D p3D, Model m3D, float *T);
};

# endif
