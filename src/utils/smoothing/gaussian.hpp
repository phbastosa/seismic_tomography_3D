# ifndef GAUSSIAN_HPP
# define GAUSSIAN_HPP

class Gaussian
{
public:

    int xdim;
    int ydim;
    int zdim;

    float stdv;           
    int samples;         

    float * volume;    

    void gaussian();
};

# endif
