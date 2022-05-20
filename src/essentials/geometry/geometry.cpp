# include "geometry.hpp"

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

void Geometry3D::setOBNS(Utils::point2D SW, Utils::point2D NW, Utils::point2D SE, float nsx, float nsy, float depth)
{
    std::vector<float> x;
    std::vector<float> y;

    this->ns = nsx * nsy;

    this->shots = new position[this->ns];

    this->shots->x = new float[this->ns];
    this->shots->y = new float[this->ns];
    this->shots->z = new float[this->ns];

    if (nsx == 1 && nsy == 1)
    {
        this->shots->x[0] = SW.x;
        this->shots->y[0] = SW.y;
        this->shots->z[0] = depth;
    }
    else if (nsx == 1)
    {
        float yi = SW.y;
        float yf = NW.y;

        y = Utils::linspace(yi, yf, nsy);

        for (int p = 0; p < nsy; p++) 
        {
            this->shots->x[p] = SW.x;   
            this->shots->y[p] = y[p];
            this->shots->z[p] = depth;    
        }
    }
    else if (nsy == 1)
    {
        float xi = SW.x;
        float xf = SE.x;

        x = Utils::linspace(xi, xf, nsx);

        for (int p = 0; p < nsx; p++) 
        {
            this->shots->x[p] = x[p];   
            this->shots->y[p] = SW.y;
            this->shots->z[p] = depth;    
        }
    }
    else
    {
        float xi = SW.x;
        float xf = SE.x;
        
        float yi = SW.y;
        float yf = NW.y;

        x = Utils::linspace(xi, xf, nsx);
        y = Utils::linspace(yi, yf, nsy);

        for (int k = 0; k < y.size(); k++)
        {
            for (int j = 0; j < x.size(); j++)
            {
                this->shots->x[j + k*x.size()] = x[j];
                this->shots->y[j + k*x.size()] = y[k];
                this->shots->z[j + k*x.size()] = depth;
            }
        }    
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setOBNR(Utils::point2D SW, Utils::point2D NW, Utils::point2D SE, float nrx, float nry, float depth)
{
    std::vector<float> x;
    std::vector<float> y;

    this->nr = nrx * nry;     

    this->nodes = new position[this->nr]; 

    this->nodes->x = new float[this->nr];
    this->nodes->y = new float[this->nr];
    this->nodes->z = new float[this->nr];

    if (nrx == 1 && nry == 1)
    {
        this->nodes->x[0] = SW.x;
        this->nodes->y[0] = SW.y;
        this->nodes->z[0] = depth;
    }
    else if (nrx == 1)
    {
        float yi = SW.y;
        float yf = NW.y;

        y = Utils::linspace(yi, yf, nry);

        for (int p = 0; p < nry; p++) 
        {
            this->nodes->x[p] = SW.x;   
            this->nodes->y[p] = y[p];
            this->nodes->z[p] = depth;    
        }
    }
    else if (nry == 1)
    {
        float xi = SW.x;
        float xf = SE.x;

        x = Utils::linspace(xi, xf, nrx);

        for (int p = 0; p < nrx; p++) 
        {
            this->nodes->x[p] = x[p];   
            this->nodes->y[p] = SW.y;
            this->nodes->z[p] = depth;    
        }
    }
    else
    {
        float xi = SW.x;
        float xf = SE.x;
        
        float yi = SW.y;
        float yf = NW.y;

        x = Utils::linspace(xi, xf, nrx);
        y = Utils::linspace(yi, yf, nry);

        for (int k = 0; k < y.size(); k++)
        {
            for (int j = 0; j < x.size(); j++)
            {
                this->nodes->x[j + k*x.size()] = x[j];
                this->nodes->y[j + k*x.size()] = y[k];
                this->nodes->z[j + k*x.size()] = depth;
            }
        }    
    }

    std::vector< float >().swap(x);
    std::vector< float >().swap(y);
}

void Geometry3D::setReciprocity()
{
    float * x = new float[this->ns];
    float * y = new float[this->ns];

    for (int p = 0; p < this->ns; p++)
    {
        x[p] = this->shots->x[p];
        y[p] = this->shots->y[p];
    }    

    delete[] this->shots->x;
    delete[] this->shots->y;

    this->shots->x = new float[this->nr];
    this->shots->y = new float[this->nr];

    for (int p = 0; p < this->nr; p++)
    {
        this->shots->x[p] = this->nodes->x[p];
        this->shots->y[p] = this->nodes->y[p];
    }    

    delete[] this->nodes->x;
    delete[] this->nodes->y;

    this->nodes->x = new float[this->ns];
    this->nodes->y = new float[this->ns];

    for (int p = 0; p < this->ns; p++)
    {
        this->nodes->x[p] = x[p];
        this->nodes->y[p] = y[p];
    }    

    int aux = this->ns; this->ns = this->nr; this->nr = aux;   

    delete[] x;
    delete[] y;
}

void Geometry3D::exportPositions(std::string shotsPath, std::string nodesPath)
{
    std::ofstream shotsFile(shotsPath, std::ios::out);        
    std::ofstream nodesFile(nodesPath, std::ios::out);

    for (int node = 0; node < this->nr; node++)        
    {   
        nodesFile <<this->nodes->x[node]<<", "<<this->nodes->y[node]<<", "<<this->nodes->z[node]<<std::endl;
    }

    for (int shot = 0; shot < this->ns; shot++)        
    {   
        shotsFile <<this->shots->x[shot]<<", "<<this->shots->y[shot]<<", "<<this->shots->z[shot]<<std::endl;
    }

    shotsFile.close();
    nodesFile.close();
}
