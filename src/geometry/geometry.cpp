# include <fstream>

# include "geometry.hpp"

void Geometry::set_topography(Coordinates &obj)
{  
    if (topography) 
    {    
        std::vector<std::string> aux_topo;

        fm.read_text_file(topography_file, aux_topo);

        for (int i = 0; i < obj.all; i++)
        {
            obj.z[i] = std::stof(aux_topo[i]);
        }

        std::vector < std::string >().swap(aux_topo);
    }
    else
    {
        for (int i = 0; i < obj.all; i++)
        {
            obj.z[i] = elevation;
        }
    }
}

void Geometry::export_positions(Coordinates &obj, std::string file)
{
    std::ofstream write(folder + file, std::ios::out);        

    for (int i = 0; i < obj.all; i++)        
    {   
        write <<obj.x[i]<<", "<<obj.y[i]<<", "<<obj.z[i]<<std::endl;    
    }

    write.close();
}

