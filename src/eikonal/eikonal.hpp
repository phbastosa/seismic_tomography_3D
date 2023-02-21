# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <string>

# include "../model/model.hpp"
# include "../model/eikonal_model/eikonal_model.hpp"

# include "../geometry/geometry.hpp"
# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

# include "../utils/file_manager/file_manager.hpp"
# include "../utils/interpolation/trilinear.hpp"

class Eikonal 
{   
private:

    File_manager fm;
    Trilinear interpolate;

protected:

    bool reciprocity;

    float * slowness;
    float * travel_time;
    float * illumination;      
    float * first_arrival;

    bool export_time_volume;  
    bool export_illumination; 
    bool export_ray_position; 
    bool export_first_arrival;

    std::string ray_folder;         
    std::string time_volume_folder; 
    std::string illumination_folder; 
    std::string first_arrival_folder;

public:

    Eikonal_model eiko_m;

    Geometry * geometry[2];

    int shot_id;
    int shots_type;
    int nodes_type;

    void info_message();
    void set_parameters(std::string file);

    void ray_tracing();
    void write_time_volume();
    void write_illumination();
    void write_first_arrival();

    virtual void solve() = 0;
    virtual void prepare_volumes() = 0;    
};

# endif
