# ifndef EIKONAL_MODEL_HPP
# define EIKONAL_MODEL_HPP

# include <string>

# include "../model.hpp"

class Eikonal_model : public Model
{
public:    
   
    float * slowness;
    float * velocity;
    float * illumination;      
    float * travel_times;

    void set_parameters(std::string file);  
};

# endif
