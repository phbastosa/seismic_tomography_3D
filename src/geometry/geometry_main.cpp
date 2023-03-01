# include "geometry.hpp"
# include "regular/regular.hpp"
# include "circular/circular.hpp"

# include "../utils/file_manager/file_manager.hpp"

int main(int argc, char **argv)
{
    Geometry * geometry[] = 
    { 
        new Regular(),
        new Circular() 
    };
    
    auto fm = File_manager();

    int shots_type = std::stoi(fm.catch_parameter("shot_geometry_type", std::string(argv[1])));
    int nodes_type = std::stoi(fm.catch_parameter("node_geometry_type", std::string(argv[1])));

    bool reciprocity = fm.str2bool(fm.catch_parameter("reciprocity", std::string(argv[1])));

    geometry[shots_type]->set_parameters(std::string(argv[1]));
    geometry[nodes_type]->set_parameters(std::string(argv[1]));

    geometry[shots_type]->export_positions(geometry[shots_type]->shots, "xyz_shots.txt");
    geometry[nodes_type]->export_positions(geometry[nodes_type]->nodes, "xyz_nodes.txt");
}
