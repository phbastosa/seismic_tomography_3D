# include <cmath>
# include "geometry.hpp"

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

void Geometry3D::setGridShots()
{
    ns = nsx * nsy;

    shots = new Position[ns];

    shots->x = new float[ns];
    shots->y = new float[ns];
    shots->z = new float[ns];

    std::vector<float> x = Utils::linspace(SW.x, SE.x, nsx);
    std::vector<float> y = Utils::linspace(SW.y, NW.y, nsy);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            shots->x[j + k*x.size()] = x[j];
            shots->y[j + k*x.size()] = y[k];
            shots->z[j + k*x.size()] = sElev;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setGridNodes()
{
    nr = nrx * nry;     

    nodes = new Position[nr]; 

    nodes->x = new float[nr];
    nodes->y = new float[nr];
    nodes->z = new float[nr];

    std::vector<float> x = Utils::linspace(SW.x, SE.x, nrx);
    std::vector<float> y = Utils::linspace(SW.y, NW.y, nry);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            nodes->x[j + k*x.size()] = x[j];
            nodes->y[j + k*x.size()] = y[k];
            nodes->z[j + k*x.size()] = rElev;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setCircularShots()
{
    std::vector<float> x, y;

    for (float radius : circles.offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {
            theta += acos(1.0f - powf(circles.ds,2.0f)/(2.0f*powf(radius,2.0f)));    
            
            x.push_back(radius*sin(theta) + circles.xc);        
            y.push_back(radius*cos(theta) + circles.yc);        

        }
    }

    ns = x.size();

    shots = new Position[ns](); 

    shots->x = new float[ns]();
    shots->y = new float[ns]();
    shots->z = new float[ns]();

    for (int i = 0; i < ns; i++)
    {
        shots->x[i] = x[i]; 
        shots->y[i] = y[i];
        shots->z[i] = sElev;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setCircularNodes()
{
    std::vector<float> x, y;

    for (float radius : circles.offsets)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {
            theta += acos(1.0f - powf(circles.ds,2.0f)/(2.0f*powf(radius,2.0f)));    
            
            x.push_back(radius*sin(theta) + circles.xc);        
            y.push_back(radius*cos(theta) + circles.yc);        

        }
    }

    nr = x.size();

    nodes = new Position[nr](); 

    nodes->x = new float[nr]();
    nodes->y = new float[nr]();
    nodes->z = new float[nr]();

    for (int i = 0; i < nr; i++)
    {
        nodes->x[i] = x[i]; 
        nodes->y[i] = y[i];
        nodes->z[i] = rElev;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setReciprocity()
{
    float * x = new float[ns];
    float * y = new float[ns];
    float * z = new float[ns];

    for (int p = 0; p < ns; p++)
    {
        x[p] = shots->x[p];
        y[p] = shots->y[p];
        z[p] = shots->z[p];
    }    

    delete[] shots->x;
    delete[] shots->y;
    delete[] shots->z;

    shots->x = new float[nr];
    shots->y = new float[nr];
    shots->z = new float[nr];

    for (int p = 0; p < nr; p++)
    {
        shots->x[p] = nodes->x[p];
        shots->y[p] = nodes->y[p];
        shots->z[p] = nodes->z[p];
    }    

    delete[] nodes->x;
    delete[] nodes->y;
    delete[] nodes->z;

    nodes->x = new float[ns];
    nodes->y = new float[ns];
    nodes->z = new float[ns];

    for (int p = 0; p < ns; p++)
    {
        nodes->x[p] = x[p];
        nodes->y[p] = y[p];
        nodes->z[p] = z[p];
    }    

    int aux = ns; ns = nr; nr = aux;   

    delete[] x;
    delete[] y;
    delete[] z;
}

void Geometry3D::exportPositions()
{
    std::ofstream shotsFile(shotsPath);        
    std::ofstream nodesFile(nodesPath);

    for (int node = 0; node < nr; node++)        
    {   
        nodesFile <<nodes->x[node]<<", "<<nodes->y[node]<<", "<<nodes->z[node]<<std::endl;
    }

    for (int shot = 0; shot < ns; shot++)        
    {   
        shotsFile <<shots->x[shot]<<", "<<shots->y[shot]<<", "<<shots->z[shot]<<std::endl;    
    }

    shotsFile.close();
    nodesFile.close();
}
