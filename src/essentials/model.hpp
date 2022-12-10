# ifndef MODEL_HPP
# define MODEL_HPP

# include <string>

class Model
{
public:    
    int nb;                  // Boundary samples
    int nx;                  // N samples in x direction
    int ny;                  // N samples in y direction
    int nz;                  // N samples in z direction   

    float dx;                // Sample spacing in x direction
    float dy;                // Sample spacing in y direction
    float dz;                // Sample spacing in z direction

    int nxx;                 // N samples in x expanded direction
    int nyy;                 // N samples in y expanded direction 
    int nzz;                 // N samples in z expanded direction 

    int nPoints;             // Total samples in model 
    int nPointsB;            // Total samples in expanded model 

    std::string vpModelPath; // Vp model location

    /* Fill model variables with other attributes dependance, invoke before declare model volumes */
    void initialize(int nx, int ny, int nz, int nb);

    /* Returns the expanded volume based on input model properties */
    float * expandVolume(float * volume, int nx, int ny, int nz, int nb);

    /* Returns the volume in original size without boundaries */
    float * reduceVolume(float * volume, int nx, int ny, int nz, int nb); 
};

# endif
