# ifndef MODEL_HPP
# define MODEL_HPP

# include "../inout/inout.hpp" 

class Model
{
private:

public:    
    int nb;                  // Boundary samples
    int nx;                  // N samples in x direction
    int ny;                  // N samples in y direction
    int nz;                  // N samples in z direction   
    int nxx;                 // N samples in x expanded direction
    int nyy;                 // N samples in y expanded direction 
    int nzz;                 // N samples in z expanded direction 

    int nPoints;             // Total samples in model 
    int nPointsB;            // Total samples in expanded model 

    float dx;                // Sample spacing in x direction
    float dy;                // Sample spacing in y direction
    float dz;                // Sample spacing in z direction

    float * vp;              // P wave velocity volume         
    float * vs;              // S wave velocity volume
    float * rho;             // Density volume
    
    std::string vpPath;
    std::string vsPath;
    std::string rhoPath;

    void initialize();

    float * readAndExpandModel(std::string path);
};

# endif
