# ifndef TRILINEAR_HPP
# define TRILINEAR_HPP

class Trilinear
{
public:

    int new_nx, new_ny, new_nz;
     
    float new_dx, new_dy, new_dz;

    float x, y, z;
    
    float x0, x1; 
    float y0, y1;
    float z0, z1; 
    
    float c000, c001, c100, c101;
    float c010, c011, c110, c111; 

    float trilinear();
};

# endif
