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

void Regular::set_parameters(std::string file)
{
    folder = fm.catch_parameter("geometry_folder", file);

    int max_coordinate_objects = 2;
    
    std::vector<std::string> coord = {"shot", "node"};

    for (int i = 0; i < max_coordinate_objects; i++)
    {
        elevation = std::stof(fm.catch_parameter(coord[i] + "_elevation", file));
        topography = fm.str2bool(fm.catch_parameter(coord[i] + "_topography", file));
        topography_file = fm.catch_parameter(coord[i] + "_topography_file", file);
    
        n_xline = std::stoi(fm.catch_parameter(coord[i] + "_n_xline", file));
        n_yline = std::stoi(fm.catch_parameter(coord[i] + "_n_yline", file));

        splitted = fm.split(fm.catch_parameter(coord[i] + "_SW", file),',');
        SW.x = std::stof(splitted[0]);
        SW.y = std::stof(splitted[1]);

        splitted = fm.split(fm.catch_parameter(coord[i] + "_NW", file),',');
        NW.x = std::stof(splitted[0]);
        NW.y = std::stof(splitted[1]);

        splitted = fm.split(fm.catch_parameter(coord[i] + "_SE", file),',');
        SE.x = std::stof(splitted[0]);
        SE.y = std::stof(splitted[1]);

        if (i == 0) build_geometry(shots);
        if (i == 1) build_geometry(nodes);
    }
}

void Regular::build_geometry(Coordinates &obj)
{
    obj.all = n_xline * n_yline;

    obj.x = new float[obj.all];
    obj.y = new float[obj.all];
    obj.z = new float[obj.all];

    std::vector<float> x_tmp = linspace(SW.x, SE.x, n_xline);
    std::vector<float> y_tmp = linspace(SW.y, NW.y, n_yline);

    for (int k = 0; k < y_tmp.size(); k++)
    {
        for (int j = 0; j < x_tmp.size(); j++)
        {
            obj.x[j + k*x_tmp.size()] = x_tmp[j];
            obj.y[j + k*y_tmp.size()] = y_tmp[k];
        }
    }    

    set_topography(obj);

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}





