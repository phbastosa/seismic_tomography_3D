# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
private:

    typedef struct
    {
        float x;
        float y;        
        
    } Point;

    Point NW, SW, SE;

    int n_xline;
    int n_yline;

    std::vector<float> linspace(float xi, float xf, int n);        

public:  

    void build_geometry();
    void set_shot_parameters(std::string file);
    void set_node_parameters(std::string file);
};

# endif
