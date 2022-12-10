# include <iostream>

# include "../../essentials/geometry.hpp"

int main(int argc, char**argv)
{
    auto G = Geometry();

    std::cout<<"\nPreparing acquisition geometry test \n"<<std::endl;

    // Point acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    G.set_SW(5000.0f, 5500.0f); // (x, y)   
    G.set_NW(5000.0f, 5500.0f); // (x, y)   
    G.set_SE(5000.0f, 5500.0f); // (x, y)   

    G.shots.n_xline = 1;
    G.shots.n_yline = 1;
    G.shots.elevation = 15.0f;

    G.setGridShots();

    // Preparing shots position

    G.set_SW(5000.0f, 4500.0f); // (x, y)   
    G.set_NW(5000.0f, 4500.0f); // (x, y)   
    G.set_SE(5000.0f, 4500.0f); // (x, y)   

    G.nodes.n_xline = 1;
    G.nodes.n_yline = 1;
    G.nodes.elevation = 500.0f;

    G.setGridNodes();    

    G.shotsPath = "outputs/pointShot.txt";
    G.nodesPath = "outputs/pointNode.txt";

    G.exportPositions();

    std::cout<<"- point acquisition test written"<<std::endl;

    // Horizontal xline acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    G.set_SW(1000.0f, 1000.0f); // (x, y)   
    G.set_NW(1000.0f, 1000.0f); // (x, y)   
    G.set_SE(9000.0f, 1000.0f); // (x, y)   

    G.shots.n_xline = 9;
    G.shots.n_yline = 1;
    G.shots.elevation = 10.0f;

    G.setGridShots();

    // Preparing receivers position

    G.set_SW(1000.0f, 9000.0f);  // (x, y)   
    G.set_NW(1000.0f, 9000.0f); // (x, y)   
    G.set_SE(9000.0f, 9000.0f);  // (x, y)   

    G.nodes.n_xline = 9;
    G.nodes.n_yline = 1;
    G.nodes.elevation = 500.0f;

    G.setGridNodes();    

    G.shotsPath = "outputs/xLineShots.txt";
    G.nodesPath = "outputs/xLineNodes.txt";

    G.exportPositions();

    std::cout<<"- xline acquisition test written"<<std::endl;

    // Vertical yline acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    G.set_SW(1000.0f, 1000.0f); // (x, y)   
    G.set_NW(1000.0f, 9000.0f); // (x, y)   
    G.set_SE(1000.0f, 1000.0f); // (x, y)   

    G.shots.n_xline = 1;
    G.shots.n_yline = 9;
    G.shots.elevation = 10.0f;

    G.setGridShots();

    // Preparing receivers position

    G.set_SW(9000.0f, 1000.0f); // (x, y)   
    G.set_NW(9000.0f, 9000.0f); // (x, y)   
    G.set_SE(9000.0f, 1000.0f); // (x, y)   

    G.nodes.n_xline = 1;
    G.nodes.n_yline = 9;
    G.nodes.elevation = 500.0f;

    G.setGridNodes();    

    G.shotsPath = "outputs/yLineShots.txt";
    G.nodesPath = "outputs/yLineNodes.txt";

    G.exportPositions();

    std::cout<<"- yline acquisition test written"<<std::endl;

    // Full carpet acquisition test -------------------------------------------------------------------------------

    // Preparing shots position

    G.set_SW(2000.0f, 2000.0f); // (x, y)   
    G.set_NW(2000.0f, 8000.0f); // (x, y)   
    G.set_SE(8000.0f, 2000.0f); // (x, y)   

    G.shots.n_xline = 7;
    G.shots.n_yline = 7;
    G.shots.elevation = 10.0f;

    G.setGridShots();

    // Preparing receivers position

    G.set_SW(2500.0f, 2500.0f);  // (x, y)   
    G.set_NW(2500.0f, 7500.0f); // (x, y)   
    G.set_SE(7500.0f, 2500.0f); // (x, y)   

    G.nodes.n_xline = 6;
    G.nodes.n_yline = 6;
    G.nodes.elevation = 500.0f;

    G.setGridNodes();    

    G.shotsPath = "outputs/carpetShots.txt";
    G.nodesPath = "outputs/carpetNodes.txt";

    G.exportPositions();

    std::cout<<"- carpet acquisition test written"<<std::endl;

    // Reciprocal carpet acquisition test -------------------------------------------------------------------------------
    
    G.setReciprocity();

    G.shotsPath = "outputs/carpetShotsReciprocity.txt";
    G.nodesPath = "outputs/carpetNodesReciprocity.txt";

    G.exportPositions();

    std::cout<<"- reciprocal carpet acquisition test written"<<std::endl;

    // Circular acquisition test --------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Preparing shots 

    G.shots.xcenter = 2500.0f;
    G.shots.ycenter = 7500.0f;
    G.shots.elevation = 10.0f;
    G.shots.circle_spacing = 25.0f;
    G.shots.offsets = {500.0f, 1000.0f, 1500.0f, 2000.0f};

    G.setCircularShots();

    // Preparing nodes

    G.nodes.xcenter = 7500.0f;
    G.nodes.ycenter = 2500.0f;
    G.nodes.elevation = 500.0f;
    G.nodes.circle_spacing = 25.0f;
    G.nodes.offsets = {500.0f, 1000.0f, 1500.0f, 2000.0f};

    G.setCircularNodes();

    G.shotsPath = "outputs/circularShots.txt";
    G.nodesPath = "outputs/circularNodes.txt";

    G.exportPositions();


    std::cout<<"- circular acquisition test written"<<std::endl;

    return 0;
}