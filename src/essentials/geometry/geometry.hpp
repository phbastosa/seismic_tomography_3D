# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

class Geometry
{
public:
    int ns;                // Total shots in simulation
    int nr;                // Total receivers in simulation 

    int nsx, nrx;          // Total geometry unit in xline
    int nsy, nry;          // Total geometry unit in crossline   

    float sElev;           // Shots elevation
    float rElev;           // Nodes elevation  

    /* Position struct to storage coordinates */
    typedef struct
    {
        float * x;         // Position x of shots or nodes
        float * y;         // Position y of shots or nodes
        float * z;         // Position z of shots or nodes
    
        int xId;
        int yId;
        int zId;

    } Position;            

    typedef struct 
    {
        float xc;
        float yc;

        float ds;

        std::vector<float> offsets;  

    } Circles;
    
    Position * shots;      // Struct to storage x,y,z coordinates of shots
    Position * nodes;      // Struct to storage x,y,z coordinates of nodes

    Circles circles;

    Utils::point2D SW;
    Utils::point2D NW;
    Utils::point2D SE;

    std::string shotsPath; 
    std::string nodesPath;

    /* Method to set shots positions giving three points in grid */
    void setGridShots();
    
    /* Method to set nodes positions giving three points in grid */
    void setGridNodes();

    /* Method to set circular shots positions */
    void setCircularShots();

    /* Method to set circular nodes positions */
    void setCircularNodes();

    /* Switch the nodes position to shots position */
    void setReciprocity();

    /* Methood to save geometry in hard drive */
    void exportPositions();
};

# endif
