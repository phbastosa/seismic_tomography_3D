# ifndef MODEL_HPP
# define MODEL_HPP

# include <string>

class Model
{
public:    
   
    int x_samples;               
    int y_samples;
    int z_samples;

    int total_samples;
    int boundary_samples;               

    float x_spacing;
    float y_spacing;
    float z_spacing;

    float * property;               

    void expand();
    void reduce(); 

    void set_parameters(std::string file);
};

# endif
