# ifndef CIRCULAR_HPP
# define CIRCULAR_HPP

# include "../geometry.hpp"

class Circular : public Geometry
{
private:

    float xcenter;               
    float ycenter;                  
    float spacing;         

    std::vector<float> offsets;  

public:  

    void build_geometry();
    void set_shot_parameters(std::string file);
    void set_node_parameters(std::string file);
};

# endif
