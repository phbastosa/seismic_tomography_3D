# include <iostream>

# include "geometry.hpp"
# include "regular/regular.hpp"
# include "circular/circular.hpp"
# include "../file_manager/file_manager.hpp"

int main(int argc, char **argv)
{
    Geometry * shots[] = 
    { 
        new Circular(), 
        new Regular() 
    };
    
    Geometry * nodes[] = 
    {
        new Circular(),
        new Regular()         
    };

    auto fm = File_manager();

    fm.parameter_file = std::string(argv[1]);
    
    int shot_type = std::stoi(fm.catch_parameter("shots_geometry_type"));
    int node_type = std::stoi(fm.catch_parameter("nodes_geometry_type"));

    bool reciprocity = fm.str2bool(fm.catch_parameter("reciprocity"));

    if (reciprocity)
    {
        shots[node_type]->set_node_parameters(fm.parameter_file);    
        nodes[shot_type]->set_shot_parameters(fm.parameter_file);

        shots[node_type]->export_positions("xyz_shots.txt");
        nodes[shot_type]->export_positions("xyz_nodes.txt");
    }
    else
    {
        shots[shot_type]->set_shot_parameters(fm.parameter_file);
        nodes[node_type]->set_node_parameters(fm.parameter_file);

        shots[shot_type]->export_positions("xyz_shots.txt");
        nodes[node_type]->export_positions("xyz_nodes.txt");
    }
}
