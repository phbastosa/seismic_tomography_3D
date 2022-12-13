# include <cmath>
# include <vector>
# include <string>
# include <fstream>

# include "geometry.hpp"

std::vector<float> Geometry::linspace(float xi, float xf, int n)
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

void Geometry::set_SW(float x, float y) {SW.x = x; SW.y = y;}

void Geometry::set_NW(float x, float y) {NW.x = x; NW.y = y;}

void Geometry::set_SE(float x, float y) {SE.x = x; SE.y = y;}

void Geometry::setGridGeometry(Position &obj)
{
    obj.all = obj.n_xline * obj.n_yline;

    obj.x = new float[obj.all];
    obj.y = new float[obj.all];
    obj.z = new float[obj.all];

    std::vector<float> x = linspace(SW.x, SE.x, obj.n_xline);
    std::vector<float> y = linspace(SW.y, NW.y, obj.n_yline);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            obj.x[j + k*x.size()] = x[j];
            obj.y[j + k*x.size()] = y[k];
            obj.z[j + k*x.size()] = obj.elevation;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setCircularGeometry(Position &obj)
{
    std::vector<float> x, y;

    for (float radius : obj.offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {            
            x.push_back(radius*sin(theta) + obj.xcenter);        
            y.push_back(radius*cos(theta) + obj.ycenter);        

            theta += acos(1.0f - powf(obj.circle_spacing,2.0f)/(2.0f*powf(radius,2.0f)));    
        }
    }

    obj.all = x.size();

    obj.x = new float[obj.all]();
    obj.y = new float[obj.all]();
    obj.z = new float[obj.all]();

    for (int i = 0; i < obj.all; i++)
    {
        obj.x[i] = x[i]; 
        obj.y[i] = y[i];
        obj.z[i] = obj.elevation;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setReciprocity()
{
    std::swap(shots.x, nodes.x);
    std::swap(shots.y, nodes.y);
    std::swap(shots.z, nodes.z);
    std::swap(shots.all, nodes.all);
}

void Geometry::exportPositions()
{
    std::ofstream shotsFile(shotsPath);        
    std::ofstream nodesFile(nodesPath);

    for (int node = 0; node < nodes.all; node++)        
    {   
        nodesFile <<nodes.x[node]<<", "<<nodes.y[node]<<", "<<nodes.z[node]<<std::endl;
    }

    for (int shot = 0; shot < shots.all; shot++)        
    {   
        shotsFile <<shots.x[shot]<<", "<<shots.y[shot]<<", "<<shots.z[shot]<<std::endl;    
    }

    shotsFile.close();
    nodesFile.close();
}
