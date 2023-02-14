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

void Regular::set_shot_parameters(std::string file)
{
    fm.parameter_file = file;

    get_shot_parameters();

    n_xline = std::stoi(fm.catch_parameter("shots_n_xline"));
    n_yline = std::stoi(fm.catch_parameter("shots_n_yline"));

    splitted = fm.split(fm.catch_parameter("shot_SW"),',');
    SW.x = std::stof(splitted[0]);
    SW.y = std::stof(splitted[1]);

    splitted = fm.split(fm.catch_parameter("shot_NW"),',');
    NW.x = std::stof(splitted[0]);
    NW.y = std::stof(splitted[1]);

    splitted = fm.split(fm.catch_parameter("shot_SE"),',');
    SE.x = std::stof(splitted[0]);
    SE.y = std::stof(splitted[1]);

    build_geometry();
}

void Regular::set_node_parameters(std::string file)
{
    fm.parameter_file = file;

    get_node_parameters();

    n_xline = std::stoi(fm.catch_parameter("nodes_n_xline"));
    n_yline = std::stoi(fm.catch_parameter("nodes_n_yline"));

    splitted = fm.split(fm.catch_parameter("node_SW"),',');
    SW.x = std::stof(splitted[0]);
    SW.y = std::stof(splitted[1]);

    splitted = fm.split(fm.catch_parameter("node_NW"),',');
    NW.x = std::stof(splitted[0]);
    NW.y = std::stof(splitted[1]);

    splitted = fm.split(fm.catch_parameter("node_SE"),',');
    SE.x = std::stof(splitted[0]);
    SE.y = std::stof(splitted[1]);

    build_geometry();
}

void Regular::build_geometry()
{
    all = n_xline * n_yline;

    x = new float[all];
    y = new float[all];
    z = new float[all];

    std::vector<float> x_tmp = linspace(SW.x, SE.x, n_xline);
    std::vector<float> y_tmp = linspace(SW.y, NW.y, n_yline);

    for (int k = 0; k < y_tmp.size(); k++)
    {
        for (int j = 0; j < x_tmp.size(); j++)
        {
            x[j + k*x_tmp.size()] = x_tmp[j];
            y[j + k*y_tmp.size()] = y_tmp[k];
        }
    }    

    set_topography();

    std::vector< float >().swap(x_tmp);
    std::vector< float >().swap(y_tmp);
}





