# include <cmath>

# include "circular.hpp"

void Circular::set_parameters(std::string file)
{
    folder = fm.catch_parameter("geometry_folder", file);

    int max_coordinate_objects = 2;

    std::vector<std::string> coord = {"shot", "node"};

    for (int i = 0; i < max_coordinate_objects; i++)
    {
        elevation = std::stof(fm.catch_parameter(coord[i] + "_elevation", file));
        topography = fm.str2bool(fm.catch_parameter(coord[i] + "_topography", file));
        topography_file = fm.catch_parameter(coord[i] + "_topography_file", file);

        xcenter = std::stof(fm.catch_parameter(coord[i] + "_x_center", file));
        ycenter = std::stof(fm.catch_parameter(coord[i] + "_y_center", file));
        spacing = std::stof(fm.catch_parameter(coord[i] + "_spacing", file));

        splitted = fm.split(fm.catch_parameter(coord[i] + "_offsets", file),',');

        for (auto offset : splitted)
            offsets.push_back(std::stof(offset));

        if (i == 0) build_geometry(shots);
        if (i == 1) build_geometry(nodes);
 
        std::vector<float>().swap(offsets);    
    }
}

void Circular::build_geometry(Coordinates &obj)
{
    std::vector<float> x_tmp, y_tmp;

    obj.all = 0;

    for (float radius : offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {            
            x_tmp.push_back(radius*sin(theta) + xcenter);        
            y_tmp.push_back(radius*cos(theta) + ycenter);        

            theta += acos(1.0f - powf(spacing, 2.0f) / (2.0f * powf(radius, 2.0f)));    

            obj.all += 1;
        }
    }

    obj.x = new float[obj.all]();
    obj.y = new float[obj.all]();
    obj.z = new float[obj.all]();

    for (int i = 0; i < obj.all; i++)
    {
        obj.x[i] = x_tmp[i]; 
        obj.y[i] = y_tmp[i];
    }

    set_topography(obj);

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}
