# include <cmath>
# include <fstream>
# include <iostream>

# include "circular.hpp"

void Circular::set_geometry(std::string parameters, std::string name)
{
    std::vector<std::string> splitted;

    elevation = std::stof(catch_parameter(name + "_elevation", parameters));
    topo_file = catch_parameter(name + "_topography_file", parameters);
    topography = str2bool(catch_parameter(name + "_topography", parameters));
    output_file = catch_parameter(name + "_output_file", parameters);

    xcenter = std::stof(catch_parameter(name + "_xcenter", parameters));
    ycenter = std::stof(catch_parameter(name + "_ycenter", parameters));
    spacing = std::stof(catch_parameter(name + "_spacing", parameters));

    splitted = split(catch_parameter(name + "_offsets", parameters),',');

    for (auto offset : splitted)
        offsets.push_back(std::stof(offset));
        
    build_geometry();
    set_topography();
    write_geometry();
 
    std::vector<float>().swap(offsets);    
}

void Circular::build_geometry()
{
    std::vector<float> x_tmp, y_tmp;

    all = 0;

    for (float radius : offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {            
            x_tmp.push_back(radius*sin(theta) + xcenter);        
            y_tmp.push_back(radius*cos(theta) + ycenter);        

            theta += acos(1.0f - powf(spacing, 2.0f) / (2.0f * powf(radius, 2.0f)));    

            all += 1;
        }
    }

    x = new float[all]();
    y = new float[all]();
    z = new float[all]();

    for (int i = 0; i < all; i++)
    {
        x[i] = x_tmp[i]; 
        y[i] = y_tmp[i];
    }

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}

void Circular::set_topography()
{
    if (topography) 
    {    
        std::vector<std::string> aux_topo;

        read_text_file(topo_file, aux_topo);

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

void Circular::write_geometry()
{
    std::ofstream file(output_file, std::ios::out);        

    if (file.is_open()) 
    {    
        for (int i = 0; i < all; i++)        
        {   
            file <<x[i]<<", "<<y[i]<<", "<<z[i]<<std::endl;    
        }
    }
    else
    {
        throw std::invalid_argument("Error: file could not be opened!");
    }

    std::cout<<"Geometry file " + output_file + " was succesfully written."<<std::endl;

    file.close();
}
