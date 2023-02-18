# include <cmath>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "eikonal.hpp"

void Eikonal::set_parameters(std::string file)
{
    model->set_parameters(file); 

    shots_type = std::stoi(fm.catch_parameter("shots_geometry_type", file)); 
    nodes_type = std::stoi(fm.catch_parameter("nodes_geometry_type", file)); 
    reciprocity = fm.str2bool(fm.catch_parameter("reciprocity", file)); 

    geometry[0] = new Regular();
    geometry[1] = new Circular();

    geometry[shots_type]->set_parameters(file);
    geometry[shots_type]->build_geometry(geometry[shots_type]->shots);
    
    geometry[nodes_type]->set_parameters(file);
    geometry[nodes_type]->build_geometry(geometry[nodes_type]->nodes);

    if (reciprocity)
    {
        std::swap(geometry[shots_type]->shots.x, geometry[nodes_type]->nodes.x); 
        std::swap(geometry[shots_type]->shots.y, geometry[nodes_type]->nodes.y); 
        std::swap(geometry[shots_type]->shots.z, geometry[nodes_type]->nodes.z); 

        std::swap(geometry[shots_type]->shots.all, geometry[nodes_type]->nodes.all);
    }

    eikonal_type = std::stoi(fm.catch_parameter("eikonalType", file));    
    
    // exportTimesVolume = str2bool(catchParameter("exportTravelTimes", parameters));
    // exportFirstArrivals = str2bool(catchParameter("exportFirstArrivals", parameters));
    // exportRayPosition = str2bool(catchParameter("exportRayPosition", parameters));
    // exportIllumination = str2bool(catchParameter("exportIllumination", parameters));

    // raysFolder = catchParameter("raysFolder", parameters);
    // eikonalFolder = catchParameter("eikonalFolder", parameters);
    // arrivalFolder = catchParameter("arrivalFolder", parameters);
    // illuminationFolder = catchParameter("illuminationFolder", parameters);

    // shots.n_xline = std::stoi(catchParameter("xShotNumber", parameters));
    // shots.n_yline = std::stoi(catchParameter("yShotNumber", parameters));

    // nodes.n_xline = std::stoi(catchParameter("xNodeNumber", parameters));
    // nodes.n_yline = std::stoi(catchParameter("yNodeNumber", parameters));

    // geometryFolder = catchParameter("geometryFolder", parameters);

    // reciprocity = str2bool(catchParameter("reciprocity", parameters));
    
    // importPositions();

    // nb = 1;
    // nx = std::stoi(catchParameter("nx", parameters));
    // ny = std::stoi(catchParameter("ny", parameters));
    // nz = std::stoi(catchParameter("nz", parameters));
    
    // initialize();

    // dx = std::stof(catchParameter("dx", parameters));
    // dy = std::stof(catchParameter("dy", parameters));
    // dz = std::stof(catchParameter("dz", parameters));
        
    // V = expand(readBinaryFloat(catchParameter("modelPath", parameters), nPoints));
    
    // if (exportIllumination) 
    //     illumination = new float[nPointsB]();
}

// void Eikonal::write_time_volume()
// {
//     if (export_time_volume)
//     {    
//         T.reduce();
        
//         fm.write_binary_float(time_volume_folder + "eikonal_nz" + std::to_string(T.z_samples) + "_nx" + std::to_string(T.x_samples) + "_ny" + std::to_string(T.y_samples) + "_shot_" + std::to_string(shot_id+1) + ".bin", T.property, T.total_samples);
//     }
// }

// void Eikonal::write_illumination()
// {
//     if (export_illumination)
//     {    
//         illumination.reduce();
        
//         fm.write_binary_float(illumination_folder + "illumination_nz" + std::to_string(illumination.z_samples) + "_nx" + std::to_string(illumination.x_samples) + "_ny" + std::to_string(illumination.y_samples) + "_shot_" + std::to_string(shot_id) + ".bin", illumination.property, illumination.total_samples);
//     }
// }

// void::Eikonal::write_first_arrival()
// {
//     if (export_first_arrival) 
//     {   
//         int nxx = T.x_samples;
//         int nyy = T.y_samples;
//         int nzz = T.z_samples;

//         int nb = T.boundary_samples;

//         float dx = T.x_spacing;
//         float dy = T.y_spacing;
//         float dz = T.z_spacing;

//         float * firstArrivals = new float[nodes[node_type]->all]();
        
//         for (int r = 0; r < nodes[node_type]->all; r++)
//         {
//             interpolation.x = nodes[node_type]->x[r];
//             interpolation.y = nodes[node_type]->y[r];
//             interpolation.z = nodes[node_type]->z[r];

//             interpolation.x0 = floorf(interpolation.x / dx) * dx;
//             interpolation.y0 = floorf(interpolation.y / dy) * dy;
//             interpolation.z0 = floorf(interpolation.z / dz) * dz;

//             interpolation.x1 = floorf(interpolation.x/dx) * dx + dx;
//             interpolation.y1 = floorf(interpolation.y/dy) * dy + dy;
//             interpolation.z1 = floorf(interpolation.z/dz) * dz + dz;

//             int id = ((int)(interpolation.z / dz) + nb) + 
//                      ((int)(interpolation.x / dx) + nb)*nzz + 
//                      ((int)(interpolation.y / dy) + nb)*nxx*nzz;

//             interpolation.c000 = T.property[id];
//             interpolation.c001 = T.property[id + 1];
//             interpolation.c100 = T.property[id + nzz]; 
//             interpolation.c101 = T.property[id + 1 + nzz]; 
//             interpolation.c010 = T.property[id + nxx*nzz]; 
//             interpolation.c011 = T.property[id + 1 + nxx*nzz]; 
//             interpolation.c110 = T.property[id + nzz + nxx*nzz]; 
//             interpolation.c111 = T.property[id + 1 + nzz + nxx*nzz];

//             firstArrivals[r] = interpolation.trilinear();        
//         }

//         writeBinaryFloat(arrivalFolder + "times_nr" + std::to_string(nodes.all) + "_shot_" + std::to_string(shotId+1) + ".bin", firstArrivals, nodes.all);

//         delete[] firstArrivals;
//     }
// }

// void Eikonal::run_ray_tracing()
// {
//     if (export_ray_position || export_illumination)
//     {
//         int nxx = T.x_samples;
//         int nyy = T.y_samples;
//         int nzz = T.z_samples;

//         int nb = T.boundary_samples;

//         float dx = T.x_spacing;
//         float dy = T.y_spacing;
//         float dz = T.z_spacing;
    
//         int sIdz = (int)(shots[shot_type]->z[shot_id] / dz) + nb;
//         int sIdx = (int)(shots[shot_type]->x[shot_id] / dx) + nb;
//         int sIdy = (int)(shots[shot_type]->y[shot_id] / dy) + nb;

//         int sId = sIdz + sIdx*nzz + sIdy*nxx*nzz;     
        
//         std::vector<float> xRay, yRay, zRay, iRay;

//         int im, jm, km, rId;

//         float rayStep = 0.2f * (dx + dy + dz) / 3.0f;

//         for (int rayId = 0; rayId < nodes[node_type]->all; rayId++)
//         {
//             float zi = nodes[node_type]->z[rayId];
//             float xi = nodes[node_type]->x[rayId];
//             float yi = nodes[node_type]->y[rayId];

//             im = (int)(zi / dz) + nb; 
//             jm = (int)(xi / dx) + nb; 
//             km = (int)(yi / dy) + nb; 

//             rId = im + jm*nzz + km*nxx*nzz;

//             if (export_illumination) illumination.property[rId] += rayStep;

//             if (export_ray_position)
//             {
//                 xRay.push_back(xi);
//                 yRay.push_back(yi);
//                 zRay.push_back(zi);
//                 iRay.push_back((float) rayId);
//             }

//             while (true)
//             {
//                 int i = (int)(zi / dz) + nb;
//                 int j = (int)(xi / dx) + nb;
//                 int k = (int)(yi / dy) + nb;

//                 float dTz = (T.property[(i+1) + j*nzz + k*nxx*nzz] - T.property[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*dz);    
//                 float dTx = (T.property[i + (j+1)*nzz + k*nxx*nzz] - T.property[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*dx);    
//                 float dTy = (T.property[i + j*nzz + (k+1)*nxx*nzz] - T.property[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*dy);

//                 float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

//                 zi -= rayStep*dTz / norm; // z ray position update   
//                 xi -= rayStep*dTx / norm; // x ray position update   
//                 yi -= rayStep*dTy / norm; // y ray position update   

//                 im = (int)(zi / dz) + nb; 
//                 jm = (int)(xi / dx) + nb; 
//                 km = (int)(yi / dy) + nb; 

//                 rId = im + jm*nzz + km*nxx*nzz;

//                 if (rId == sId)
//                 {
//                     if (exportRayPosition)
//                     {
//                         xRay.push_back(shots.x[shotId]);
//                         yRay.push_back(shots.y[shotId]);
//                         zRay.push_back(shots.z[shotId]);
//                         iRay.push_back((float) rayId);
//                     }
                
//                     if (exportIllumination)
//                     {
//                         float finalDist = sqrtf(powf(shots.x[shotId] - xi, 2.0f) + powf(shots.y[shotId] - yi, 2.0f) + powf(shots.z[shotId] - zi, 2.0f));        
//                         illumination[rId] += finalDist;                
//                     }            

//                     break;
//                 }

//                 if (exportIllumination) illumination[rId] += rayStep;

//                 if (exportRayPosition)
//                 {
//                     xRay.push_back(xi);
//                     yRay.push_back(yi);
//                     zRay.push_back(zi);
//                     iRay.push_back((float) rayId);
//                 }
//             }
//         }

//         if (exportRayPosition)
//         {
//             std::ofstream fout;

//             size_t allRayPoints = iRay.size();

//             std::string xRayPath = raysFolder + "xRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
//             std::string yRayPath = raysFolder + "yRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
//             std::string zRayPath = raysFolder + "zRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";
//             std::string iRayPath = raysFolder + "iRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shotId+1) + ".bin";

//             fout.open(xRayPath, std::ios::out);
//             fout.write((char*)&xRay[0], xRay.size() * sizeof(xRay));
//             std::cout<<"Ray x position file was written with name " + xRayPath<<std::endl;  
//             fout.close();

//             fout.open(yRayPath, std::ios::out);
//             fout.write((char*)&yRay[0], yRay.size() * sizeof(yRay));
//             std::cout<<"Ray y position file was written with name " + yRayPath<<std::endl;  
//             fout.close();

//             fout.open(zRayPath, std::ios::out);
//             fout.write((char*)&zRay[0], zRay.size() * sizeof(zRay));
//             std::cout<<"Ray z position file was written with name " + zRayPath<<std::endl;  
//             fout.close();

//             fout.open(iRayPath, std::ios::out);
//             fout.write((char*)&iRay[0], iRay.size() * sizeof(iRay));
//             std::cout<<"Ray node index file was written with name " + iRayPath<<std::endl;  
//             fout.close();
//         }

//         std::vector<float>().swap(iRay);
//         std::vector<float>().swap(xRay);
//         std::vector<float>().swap(yRay);
//         std::vector<float>().swap(zRay);
//     }
// }

