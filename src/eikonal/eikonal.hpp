# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <string>

# include "../geometry/geometry.hpp"
# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"
// # include "../geometry/streamer/streamer.hpp"
// # include "../geometry/crosswell/crosswell.hpp"

# include "../utils/file_manager/file_manager.hpp"
# include "../utils/interpolation/trilinear.hpp"

class Eikonal 
{   
protected:

    int nx, ny, nz;
    int nxx, nyy, nzz;
    int nPoints, nPointsB;
    int total_nodes;

    float dx, dy, dz;

    float * slowness;
    float * travel_time;    
    float * first_arrival;

    bool reciprocity;
    bool export_time_volume;  
    bool export_first_arrival;

    std::string time_volume_folder;          
    std::string first_arrival_folder;

    Geometry * shots;
    Geometry * nodes; 

    virtual void expand_model() = 0;
    virtual void reduce_model() = 0;

public:

    std::string parameters;

    int shot_id;
    int total_shots;

    virtual void set_parameters() = 0;
    virtual void prepare_volumes() = 0;

    virtual void info_message() = 0;
    virtual void eikonal_equation() = 0;
    virtual void write_time_volume() = 0;
    virtual void write_first_arrival() = 0;

    virtual void destroy_volumes() = 0;    
};

# endif
