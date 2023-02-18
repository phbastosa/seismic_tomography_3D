# include "eikonal_model.hpp"

void Eikonal_model::set_parameters(std::string file)
{
    auto fm = File_manager();

    x_samples = std::stoi(fm.catch_parameter("x_samples", file));
    y_samples = std::stoi(fm.catch_parameter("y_samples", file));
    z_samples = std::stoi(fm.catch_parameter("z_samples", file));

    boundary_samples = 1;        

    x_samples_b = x_samples + 2 * boundary_samples;
    y_samples_b = y_samples + 2 * boundary_samples;
    z_samples_b = z_samples + 2 * boundary_samples;

    total_samples = x_samples * y_samples * z_samples;
    total_samples_b = x_samples_b * y_samples_b * z_samples_b;

    x_spacing = std::stof(fm.catch_parameter("x_spacing", file));    
    y_spacing = std::stof(fm.catch_parameter("y_spacing", file));    
    z_spacing = std::stof(fm.catch_parameter("z_spacing", file));    

    velocity = new float[total_samples]();
    illumination = new float[total_samples]();
    travel_times = new float[total_samples_b]();

    std::string vp_file = fm.catch_parameter("vp_file", file);

    fm.read_binary_float(vp_file, velocity, total_samples);

    slowness = expand(velocity); 

    for (int index = 0; index < total_samples_b; index++)
        slowness[index] = 1.0f / slowness[index];
}


