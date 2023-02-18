# ifndef MODEL_HPP
# define MODEL_HPP

# include "../utils/smoothing/gaussian.hpp"
# include "../utils/interpolation/trilinear.hpp"

class Model
{
private:

    Gaussian smoother;
    Trilinear interpolate;

public:    
   
    int x_samples;               
    int y_samples;
    int z_samples;

    int x_samples_b;
    int y_samples_b;
    int z_samples_b;

    int total_samples;
    int total_samples_b;
    int boundary_samples;               

    float x_spacing;
    float y_spacing;
    float z_spacing;

    int new_x_samples;
    int new_y_samples;
    int new_z_samples;

    float new_x_spacing;
    float new_y_spacing;
    float new_z_spacing; 

    int smoothing_samples;
    float standard_deviation;

    float * expand(float * volume);
    float * reduce(float * volume); 
    float * smooth(float * volume);
    float * resize(float * volume);
};

# endif
