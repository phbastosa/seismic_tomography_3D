# ifndef MODEL_HPP
# define MODEL_HPP

# include <string>

# include "../../utils/file_manager/file_manager.hpp"

class Model
{
protected:

    int smoothing_samples;
    float standard_deviation;

    int new_x_samples;
    int new_y_samples;
    int new_z_samples;

    float new_x_spacing;
    float new_y_spacing;
    float new_z_spacing;    

    File_manager fm;

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

    float * expand(float * volume);
    float * reduce(float * volume); 
    float * smooth(float * volume);
    float * resize(float * volume);

    virtual void set_parameters(std::string file) = 0;
};

# endif
