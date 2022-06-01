# include <cmath>
# include "geometry.hpp"

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

void Geometry::setGridShots()
{
    shots.n = shots.nx * shots.ny;

    shots.x = new float[shots.n];
    shots.y = new float[shots.n];
    shots.z = new float[shots.n];

    std::vector<float> x = Utils::linspace(SW.x, SE.x, shots.nx);
    std::vector<float> y = Utils::linspace(SW.y, NW.y, shots.ny);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            shots.x[j + k*x.size()] = x[j];
            shots.y[j + k*x.size()] = y[k];
            shots.z[j + k*x.size()] = shots.elevation;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setGridNodes()
{
    nodes.n = nodes.nx * nodes.ny;     

    nodes.x = new float[nodes.n];
    nodes.y = new float[nodes.n];
    nodes.z = new float[nodes.n];

    std::vector<float> x = Utils::linspace(SW.x, SE.x, nodes.nx);
    std::vector<float> y = Utils::linspace(SW.y, NW.y, nodes.ny);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            nodes.x[j + k*x.size()] = x[j];
            nodes.y[j + k*x.size()] = y[k];
            nodes.z[j + k*x.size()] = nodes.elevation;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setCircularShots()
{
    std::vector<float> x, y;

    for (float radius : shots.offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {            
            x.push_back(radius*sin(theta) + shots.xc);        
            y.push_back(radius*cos(theta) + shots.yc);        

            theta += acos(1.0f - powf(shots.ds,2.0f)/(2.0f*powf(radius,2.0f)));    
        }
    }

    shots.n = x.size();

    shots.x = new float[shots.n]();
    shots.y = new float[shots.n]();
    shots.z = new float[shots.n]();

    for (int i = 0; i < shots.n; i++)
    {
        shots.x[i] = x[i]; 
        shots.y[i] = y[i];
        shots.z[i] = shots.elevation;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setCircularNodes()
{
    std::vector<float> x, y;

    for (float radius : nodes.offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {
            x.push_back(radius*sin(theta) + nodes.xc);        
            y.push_back(radius*cos(theta) + nodes.yc);        

            theta += acos(1.0f - powf(nodes.ds,2.0f)/(2.0f*powf(radius,2.0f)));    
        }
    }

    nodes.n = x.size();

    nodes.x = new float[nodes.n]();
    nodes.y = new float[nodes.n]();
    nodes.z = new float[nodes.n]();

    for (int i = 0; i < nodes.n; i++)
    {
        nodes.x[i] = x[i]; 
        nodes.y[i] = y[i];
        nodes.z[i] = nodes.elevation;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry::setReciprocity()
{
    float * x = new float[shots.n];
    float * y = new float[shots.n];
    float * z = new float[shots.n];

    for (int p = 0; p < shots.n; p++)
    {
        x[p] = shots.x[p];
        y[p] = shots.y[p];
        z[p] = shots.z[p];
    }    

    delete[] shots.x;
    delete[] shots.y;
    delete[] shots.z;

    shots.x = new float[nodes.n];
    shots.y = new float[nodes.n];
    shots.z = new float[nodes.n];

    for (int p = 0; p < nodes.n; p++)
    {
        shots.x[p] = nodes.x[p];
        shots.y[p] = nodes.y[p];
        shots.z[p] = nodes.z[p];
    }    

    delete[] nodes.x;
    delete[] nodes.y;
    delete[] nodes.z;

    nodes.x = new float[shots.n];
    nodes.y = new float[shots.n];
    nodes.z = new float[shots.n];

    for (int p = 0; p < shots.n; p++)
    {
        nodes.x[p] = x[p];
        nodes.y[p] = y[p];
        nodes.z[p] = z[p];
    }    

    int aux = shots.n; shots.n = nodes.n; nodes.n = aux;   

    delete[] x;
    delete[] y;
    delete[] z;
}

void Geometry::exportPositions()
{
    std::ofstream shotsFile(shotsPath);        
    std::ofstream nodesFile(nodesPath);

    for (int node = 0; node < nodes.n; node++)        
    {   
        nodesFile <<nodes.x[node]<<", "<<nodes.y[node]<<", "<<nodes.z[node]<<std::endl;
    }

    for (int shot = 0; shot < shots.n; shot++)        
    {   
        shotsFile <<shots.x[shot]<<", "<<shots.y[shot]<<", "<<shots.z[shot]<<std::endl;    
    }

    shotsFile.close();
    nodesFile.close();
}
