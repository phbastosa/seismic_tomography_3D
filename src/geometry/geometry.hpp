# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include <string>

# include "../utils/file_manager/file_manager.hpp"

class Geometry
{
protected:

    bool topography;
    float elevation;
    std::string topo_file;
    std::string output_file;

    virtual void set_topography() = 0;
    virtual void build_geometry() = 0;
    virtual void write_geometry() = 0;

public:  

    int all;

    float * x;
    float * y;
    float * z;

    virtual void set_geometry(std::string parameters, std::string name) = 0;
};

# endif
