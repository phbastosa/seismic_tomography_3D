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

    void set_topography();
    void build_geometry();
    void write_geometry();

public:  

    void set_geometry(std::string parameters, std::string name);
};

# endif
