# include <fstream>
# include <iostream>

# include "regular.hpp"

std::vector<float> Regular::linspace(float xi, float xf, int n)
{
    std::vector<float> linspaced;
    
    if (n == 0) return linspaced;
    if (n == 1)
    {
        linspaced.push_back(xi);
        return linspaced;
    } 

    linspaced.reserve(n);

    float delta = (xf - xi) / (n - 1);

    for (int i = 0; i < n; i++)
    {
        linspaced.emplace_back(xi + (float)(delta*i));
    }

    return linspaced;
}

void Regular::set_geometry(std::string parameters, std::string name)
{
    std::vector<std::string> splitted;

    elevation = std::stof(catch_parameter(name + "_elevation", parameters));
    topography = str2bool(catch_parameter(name + "_topography", parameters));
    topo_file = catch_parameter(name + "_topography_file", parameters);
    output_file = catch_parameter(name + "_output_file", parameters);

    n_xline = std::stoi(catch_parameter(name + "_n_xline", parameters));
    n_yline = std::stoi(catch_parameter(name + "_n_yline", parameters));

    splitted = split(catch_parameter(name + "_SW", parameters),',');
    SW.x = std::stof(splitted[0]);
    SW.y = std::stof(splitted[1]);

    splitted = split(catch_parameter(name + "_NW", parameters),',');
    NW.x = std::stof(splitted[0]);
    NW.y = std::stof(splitted[1]);

    splitted = split(catch_parameter(name + "_SE", parameters),',');
    SE.x = std::stof(splitted[0]);
    SE.y = std::stof(splitted[1]);
        
    build_geometry();
    set_topography();
    write_geometry();
}

void Regular::set_topography()
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

void Regular::build_geometry()
{
    all = n_xline * n_yline;

    x = new float[all]();
    y = new float[all]();
    z = new float[all]();

    std::vector<float> x_tmp = linspace(SW.x, SE.x, n_xline);
    std::vector<float> y_tmp = linspace(SW.y, NW.y, n_yline);

    for (int k = 0; k < y_tmp.size(); k++)
    {
        for (int j = 0; j < x_tmp.size(); j++)
        {
            x[k + j*n_yline] = x_tmp[j];
            y[k + j*n_yline] = y_tmp[k];
        }
    }    

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}

void Regular::write_geometry()
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

