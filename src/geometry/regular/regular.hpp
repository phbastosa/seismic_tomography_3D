# ifndef REGULAR_HPP
# define REGULAR_HPP

# include <vector>

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

    void set_topography();
    void build_geometry();
    void write_geometry();

public:  

    void set_geometry(std::string parameters, std::string name);
};

# endif