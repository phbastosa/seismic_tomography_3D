# include <fstream>

# include "geometry.hpp"

void Geometry::get_shot_parameters()
{
    elevation = std::stof(fm.catch_parameter("shots_elevation"));
    topography = fm.str2bool(fm.catch_parameter("shots_topography"));
    topography_file = fm.catch_parameter("shots_topography_file");
    geometry_folder = fm.catch_parameter("geometry_folder");
}

void Geometry::get_node_parameters()
{
    elevation = std::stof(fm.catch_parameter("nodes_elevation"));
    topography = fm.str2bool(fm.catch_parameter("nodes_topography"));
    topography_file = fm.catch_parameter("nodes_topography_file");
    geometry_folder = fm.catch_parameter("geometry_folder");
}

void Geometry::set_topography()
{  
    if (topography) 
    {    
        std::vector<std::string> aux_topo;

        fm.read_text_file(topography_file, aux_topo);

        for (int i = 0; i < all; i++)
        {
            z[i] = std::stof(aux_topo[i]);
        }

        std::vector < std::string >().swap(aux_topo);
    }
    else
    {
        for (int i = 0; i < all; i++)
        {
            z[i] = elevation;
        }
    }
}

void Geometry::import_positions(std::string file)
{
    std::vector<std::string> xyz;

    fm.read_text_file(geometry_folder + file, xyz);
    
    all = xyz.size();

    for (int i = 0; i < all; i++)
    {
        splitted = fm.split(xyz[i], ',');
        
        x[i] = stof(splitted[0]);
        y[i] = stof(splitted[1]);
        z[i] = stof(splitted[2]);
    }
}

void Geometry::export_positions(std::string file)
{
    std::string txt = geometry_folder + file;

    std::ofstream write(txt, std::ios::out);        

    for (int i = 0; i < all; i++)        
    {   
        write <<x[i]<<", "<<y[i]<<", "<<z[i]<<std::endl;    
    }

    write.close();
}

