# ifndef MODEL_HPP
# define MODEL_HPP

# include <string>

class Model
{
public:    
   
    int nb;               // Boundary samples
    int nx;               // N samples in x direction
    int ny;               // N samples in y direction
    int nz;               // N samples in z direction   

    int nPoints;          // Total samples in model 

    float dx;             // Sample spacing in x direction
    float dy;             // Sample spacing in y direction
    float dz;             // Sample spacing in z direction

    float * value;        // Value of property in model       

    void expand();
    void reduce(); 
    void smooth(float stdv, int samples);
    void resize(float new_dz, float new_dx, float new_dy);

    void set_parameters(std::string file);
};

# endif
