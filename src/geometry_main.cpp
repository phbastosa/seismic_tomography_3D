# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"

int main(int argc, char **argv)
{
    Geometry * shots[] = 
    { 
        new Regular(),
        new Circular(),
        // Streamer(),
        // Crosswell() 
    };

    Geometry * nodes[] = 
    { 
        new Regular(),
        new Circular(),
        // Streamer(),
        // Crosswell() 
    };

    int shots_type = std::stoi(catch_parameter("shots_geometry_type", std::string(argv[1])));
    int nodes_type = std::stoi(catch_parameter("nodes_geometry_type", std::string(argv[1])));

    shots[shots_type]->set_geometry(std::string(argv[1]), "shots");
    nodes[nodes_type]->set_geometry(std::string(argv[1]), "nodes");
}
