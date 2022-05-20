# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../inout/inout.hpp"
# include "../utils/utils.hpp"

class Geometry3D
{
private:


public:
    int ns;                // Total shots in simulation
    int nr;                // Total receivers in simulation 

    /* Position struct to storage coordinates */
    typedef struct
    {
        float * x;         // Position x of shots or nodes
        float * y;         // Position y of shots or nodes
        float * z;         // Position z of shots or nodes
    
    } position;            

    position * shots;      // Struct to storage x,y,z coordinates of shots
    position * nodes;      // Struct to storage x,y,z coordinates of nodes

    /* Method to set shots positions giving three points in grid */
    void setOBNS(Utils::point2D SW, Utils::point2D NW, Utils::point2D SE, float nsx, float nsy, float depth);
    
    /* Method to set shots positions giving three points in grid */
    void setOBNR(Utils::point2D SW, Utils::point2D NW, Utils::point2D SE, float nsx, float nsy, float depth);

    /* Switch the nodes position to shots position */
    void setReciprocity();

    /* Methood to save geometry in hard drive */
    void exportPositions(std::string shotsPath, std::string nodesPath);
};

# endif
