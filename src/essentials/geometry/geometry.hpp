# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

class Geometry
{
public:  

    /* Position struct to storage coordinates */
    typedef struct
    {
        float * x;         // Position x of shots or nodes
        float * y;         // Position y of shots or nodes
        float * z;         // Position z of shots or nodes
    
        int idx;
        int idy;
        int idz;

        int n;             // Total elements in simulation
        int nx;            // Total geometry unit in xline
        int ny;            // Total geometry unit in crossline

        float xc;
        float yc;

        float ds;

        float elevation;

        std::vector<float> offsets;  

    } Position;            
    
    Position shots;      // Struct to storage x,y,z coordinates of shots
    Position nodes;      // Struct to storage x,y,z coordinates of nodes

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
