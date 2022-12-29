# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include <cmath>
# include <vector>
# include <string>
# include <fstream>

typedef struct                   // Position struct to storage coordinates    
{
    float * x;                   // Position x of geometry element
    float * y;                   // Position y of geometry element
    float * z;                   // Position z of geometry element

    int n;

} Position;

typedef struct                   // Struct to storage a simple 2D point
{
    float x;                     // x plane position
    float y;                     // y plane position

} Point;            

/* Function to calculate a linear spaced position */
std::vector<float> linspace(float xi, float xf, int n)
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

/* Method to set shots positions giving three points in grid */
Position setGridGeometry(Point SE, Point SW, Point NW, int n_xline, int n_yline, float elev)
{
    Position obj;

    obj.n = n_xline * n_yline;

    obj.x = new float[obj.n];
    obj.y = new float[obj.n];
    obj.z = new float[obj.n];

    std::vector<float> x = linspace(SW.x, SE.x, n_xline);
    std::vector<float> y = linspace(SW.y, NW.y, n_yline);

    for (int k = 0; k < y.size(); k++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            obj.x[j + k*x.size()] = x[j];
            obj.y[j + k*x.size()] = y[k];
            obj.z[j + k*x.size()] = elev;
        }
    }    

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);

    return obj;
}

/* Method to set nodes positions giving three points in grid */
Position setCircularGeometry(std::vector<float> circleRadius, float circleSpacing, float xCenter, float yCenter, float elev)
{
    Position obj;

    std::vector<float> x, y;

    for (float radius : circleRadius)
    {
        float theta = 0.0f;

        while (theta < 2.0f * 4.0f*atan(1.0f))
        {            
            x.push_back(radius*sin(theta) + xCenter);        
            y.push_back(radius*cos(theta) + yCenter);        

            theta += acos(1.0f - powf(circleSpacing, 2.0f)/(2.0f*powf(radius, 2.0f)));    
        }
    }

    obj.n = x.size();

    obj.x = new float[obj.n]();
    obj.y = new float[obj.n]();
    obj.z = new float[obj.n]();

    for (int i = 0; i < obj.n; i++)
    {
        obj.x[i] = x[i]; 
        obj.y[i] = y[i];
        obj.z[i] = elev;
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);

    return obj;
}

/* Switch the nodes position to shots position */
void setReciprocity(Position shots, Position nodes)
{
    std::swap(shots.x, nodes.x);
    std::swap(shots.y, nodes.y);
    std::swap(shots.z, nodes.z);
    std::swap(shots.n, nodes.n);
}

/* Method to save geometry in hard drive */
void exportPositions(std::string folder, Position shots, Position nodes)
{
    std::ofstream shotsFile(folder + "shots.txt");        
    std::ofstream nodesFile(folder + "nodes.txt");

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

# endif
