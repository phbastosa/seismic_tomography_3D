# include <iostream>

# include "../../essentials/geometry.hpp"

int main(int argc, char**argv)
{
    std::cout<<"Preparing acquisition geometry test \n"<<std::endl;

    // Point acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    Point SW, NW, SE;

    SE.x = 5000.0f; SE.y = 5500.0f; // (x, y)   NW   
    NW.x = 5000.0f; NW.y = 5500.0f; // (x, y)   | 
    SW.x = 5000.0f; SW.y = 5500.0f; // (x, y)   SW -- SE

    auto shots = setGridGeometry(SE, SW, NW, 1, 1, 15.0f);

    // Preparing nodes position

    SW.x = 5000.0f; SW.y = 4500.0f; // (x, y)   NW
    NW.x = 5000.0f; NW.y = 4500.0f; // (x, y)   |
    SE.x = 5000.0f; SE.y = 4500.0f; // (x, y)   SW -- SE

   auto nodes = setGridGeometry(SE, SW, NW, 1, 1, 25.0f);    

    exportPositions("outputs/point_", shots, nodes);

    std::cout<<"- point acquisition test written"<<std::endl;

    // Horizontal xline acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    SE.x = 1000.0f; SE.y = 1000.0f; // (x, y)   NW   
    NW.x = 1000.0f; NW.y = 1000.0f; // (x, y)   | 
    SW.x = 9000.0f; SW.y = 1000.0f; // (x, y)   SW -- SE

    shots = setGridGeometry(SE, SW, NW, 9, 1, 10.0f);

    // Preparing receivers position

    SE.x = 1000.0f; SE.y = 9000.0f; // (x, y)   NW   
    NW.x = 1000.0f; NW.y = 9000.0f; // (x, y)   | 
    SW.x = 9000.0f; SW.y = 9000.0f; // (x, y)   SW -- SE

    nodes = setGridGeometry(SE, SW, NW, 9, 1, 500.0f);

    exportPositions("outputs/xLine_", shots, nodes);

    std::cout<<"- xline acquisition test written"<<std::endl;

    // Vertical yline acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    SE.x = 1000.0f; SE.y = 1000.0f; // (x, y)   NW   
    NW.x = 1000.0f; NW.y = 9000.0f; // (x, y)   | 
    SW.x = 1000.0f; SW.y = 1000.0f; // (x, y)   SW -- SE

    shots = setGridGeometry(SE, SW, NW, 1, 9, 10.0f);

    // Preparing receivers position

    SE.x = 9000.0f; SE.y = 1000.0f; // (x, y)   NW   
    NW.x = 9000.0f; NW.y = 9000.0f; // (x, y)   | 
    SW.x = 9000.0f; SW.y = 1000.0f; // (x, y)   SW -- SE

    nodes = setGridGeometry(SE, SW, NW, 1, 9, 500.0f);

    exportPositions("outputs/yLine_", shots, nodes);

    std::cout<<"- yline acquisition test written"<<std::endl;

    // Full carpet acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    SE.x = 2000.0f; SE.y = 2000.0f; // (x, y)   NW   
    NW.x = 2000.0f; NW.y = 8000.0f; // (x, y)   | 
    SW.x = 8000.0f; SW.y = 2000.0f; // (x, y)   SW -- SE

    shots = setGridGeometry(SE, SW, NW, 7, 7, 10.0f);

    // Preparing receivers position

    SE.x = 2500.0f; SE.y = 2500.0f; // (x, y)   NW   
    NW.x = 2500.0f; NW.y = 7500.0f; // (x, y)   | 
    SW.x = 7500.0f; SW.y = 2500.0f; // (x, y)   SW -- SE

    nodes = setGridGeometry(SE, SW, NW, 6, 6, 500.0f);

    exportPositions("outputs/carpet_", shots, nodes);

    std::cout<<"- carpet acquisition test written"<<std::endl;

    // Reciprocal carpet acquisition test -------------------------------------------------------------------------------
    
    setReciprocity(shots, nodes);

    exportPositions("outputs/carpetReciprocity_",shots, nodes);

    std::cout<<"- reciprocal carpet acquisition test written"<<std::endl;

    // Circular acquisition test --------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Preparing shots 

    float spacing = 25.0f;
    float xcenter = 2500.0f;
    float ycenter = 7500.0f;
    float elevation = 10.0f;
    std::vector<float> offsets = {500.0f, 1000.0f, 1500.0f, 2000.0f};

    shots = setCircularGeometry(offsets, spacing, xcenter, ycenter, elevation);

    // Preparing nodes

    spacing = 25.0f;
    xcenter = 7500.0f;
    ycenter = 2500.0f;
    elevation = 500.0f;
    offsets = {500.0f, 1000.0f, 1500.0f, 2000.0f};

    nodes = setCircularGeometry(offsets, spacing, xcenter, ycenter, elevation);

    exportPositions("outputs/circular_", shots, nodes);

    std::cout<<"- circular acquisition test written"<<std::endl;

    return 0;
}
