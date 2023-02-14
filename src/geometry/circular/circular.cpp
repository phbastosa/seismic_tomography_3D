# include <cmath>

# include "circular.hpp"

void Circular::set_shot_parameters(std::string file)
{
    fm.parameter_file = file;

    get_shot_parameters();

    xcenter = std::stof(fm.catch_parameter("shots_xcenter"));
    ycenter = std::stof(fm.catch_parameter("shots_ycenter"));
    spacing = std::stof(fm.catch_parameter("shots_spacing"));

    splitted = fm.split(fm.catch_parameter("shots_offsets"),',');

    for (auto offset : splitted)
        offsets.push_back(std::stof(offset));

    build_geometry();    
}

void Circular::set_node_parameters(std::string file)
{
    fm.parameter_file = file;

    get_node_parameters();

    xcenter = std::stof(fm.catch_parameter("shots_xcenter"));
    ycenter = std::stof(fm.catch_parameter("shots_ycenter"));
    spacing = std::stof(fm.catch_parameter("shots_spacing"));

    splitted = fm.split(fm.catch_parameter("shots_offsets"),',');

    for (auto offset : splitted)
        offsets.push_back(std::stof(offset));

    build_geometry();
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

    set_topography();

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}
