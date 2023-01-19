# include <iostream>

# include "../../src/essentials/utils.hpp"
# include "../../src/essentials/geometry.hpp"

int main(int argc, char **argv)
{
    auto utils = Utils();    
    auto geometry = Geometry();    
    
    std::string parameters = std::string(argv[1]);

    geometry.geometryFolder = utils.catchParameter("geometryFolder", parameters);

    std::vector<std::string> splitted;

    geometry.shots.elevation = std::stof(utils.catchParameter("shotsElevation", parameters));
    geometry.nodes.elevation = std::stof(utils.catchParameter("nodesElevation", parameters));

    geometry.shots.n_xline = std::stoi(utils.catchParameter("xShotNumber", parameters));
    geometry.shots.n_yline = std::stoi(utils.catchParameter("yShotNumber", parameters));
    
    splitted = utils.split(utils.catchParameter("shotSW", parameters),',');
    geometry.set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("shotNW", parameters),',');
    geometry.set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("shotSE", parameters),',');
    geometry.set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    geometry.setGridGeometry(geometry.shots);

    geometry.nodes.n_xline = std::stoi(utils.catchParameter("xNodeNumber", parameters));
    geometry.nodes.n_yline = std::stoi(utils.catchParameter("yNodeNumber", parameters));
    
    splitted = utils.split(utils.catchParameter("nodeSW", parameters),',');
    geometry.set_SW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("nodeNW", parameters),',');
    geometry.set_NW(std::stof(splitted[0]), std::stof(splitted[1]));

    splitted = utils.split(utils.catchParameter("nodeSE", parameters),',');
    geometry.set_SE(std::stof(splitted[0]), std::stof(splitted[1]));

    geometry.setGridGeometry(geometry.nodes);

    geometry.shotsTopography = utils.str2bool(utils.catchParameter("shotsTopography", parameters));    
    geometry.nodesTopography = utils.str2bool(utils.catchParameter("nodesTopography", parameters));    

    if (geometry.shotsTopography) 
    {
        geometry.shotsTopographyPath = utils.catchParameter("shotsTopographyPath", parameters);
        geometry.shots.z = utils.readBinaryFloat(geometry.geometryFolder + geometry.shotsTopographyPath, geometry.shots.all);
    }

    if (geometry.nodesTopography)
    {
        geometry.nodesTopographyPath = utils.catchParameter("nodesTopographyPath", parameters);
        geometry.nodes.z = utils.readBinaryFloat(geometry.geometryFolder + geometry.nodesTopographyPath, geometry.nodes.all);
    }

    utils.writeBinaryFloat(geometry.geometryFolder + "x_shots_" + std::to_string(geometry.shots.all) + "_positions.bin", geometry.shots.x, geometry.shots.all);
    utils.writeBinaryFloat(geometry.geometryFolder + "y_shots_" + std::to_string(geometry.shots.all) + "_positions.bin", geometry.shots.y, geometry.shots.all);
    utils.writeBinaryFloat(geometry.geometryFolder + "z_shots_" + std::to_string(geometry.shots.all) + "_positions.bin", geometry.shots.z, geometry.shots.all);

    utils.writeBinaryFloat(geometry.geometryFolder + "x_nodes_" + std::to_string(geometry.nodes.all) + "_positions.bin", geometry.nodes.x, geometry.nodes.all);
    utils.writeBinaryFloat(geometry.geometryFolder + "y_nodes_" + std::to_string(geometry.nodes.all) + "_positions.bin", geometry.nodes.y, geometry.nodes.all);
    utils.writeBinaryFloat(geometry.geometryFolder + "z_nodes_" + std::to_string(geometry.nodes.all) + "_positions.bin", geometry.nodes.z, geometry.nodes.all);

    geometry.exportPositions();
}