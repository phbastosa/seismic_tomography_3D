# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/file_manager/file_manager.hpp"

class Geometry
{
private:

protected:

    File_manager fm;

    typedef struct 
    {
        float * x;
        float * y;
        float * z;

        int all;

    } Coordinates;

    std::vector<std::string> splitted;

    bool topography;
    float elevation;
    std::string topography_file;
    void set_topography(Coordinates &obj);

    std::string folder;

public:  
 
    Coordinates shots;
    Coordinates nodes;

    virtual void build_geometry(Coordinates &obj) = 0;
    virtual void set_parameters(std::string file) = 0;

    void export_positions(Coordinates &obj, std::string file);
};

# endif



