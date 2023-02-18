# ifndef TRILINEAR_HPP
# define TRILINEAR_HPP

class Trilinear
{
public:

    float x, y, z;
    
    float x0, x1; 
    float y0, y1;
    float z0, z1; 
    
    float c000, c001, c100, c101;
    float c010, c011, c110, c111; 

    float trilinear();
};

# endif
