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

    int eikonal_type;

    bool export_time_volume;  
    bool export_first_arrival;
    bool export_illumination; 
    bool export_ray_position; 

    std::string ray_folder;         
    std::string time_volume_folder; 
    std::string first_arrival_Folder;
    std::string illumination_folder; 
    
    File_manager fm;

protected:

    Model * model = new Eikonal_model();
    
    Geometry * geometry[2];

    int shots_type;
    int nodes_type;
    bool reciprocity;

    // Model illumination;

    // File_manager fm;

    // int shot_id;
    // int eikonal_type;

    // int shot_type;
    // int node_type;

    // bool reciprocity;

public:

    virtual void run_solver() = 0;    

    void write_time_volume();
    void write_illumination();
    void write_first_arrival();

    void run_ray_tracing();

    void set_parameters(std::string file);
};

# endif
