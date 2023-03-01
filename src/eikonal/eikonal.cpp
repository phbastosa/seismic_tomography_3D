# include <cmath>
# include <iomanip>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "eikonal.hpp"

void Eikonal::info_message()
{
    std::cout << std::fixed << std::setprecision(1);

    system("clear");
    std::cout<<"3D eikonal equation solver\n\n";
    
    std::cout<<"Total x model length = "<<(eiko_m.x_samples-1)*eiko_m.x_spacing<<" m\n";
    std::cout<<"Total Y model length = "<<(eiko_m.y_samples-1)*eiko_m.y_spacing<<" m\n";
    std::cout<<"Total Z model length = "<<(eiko_m.z_samples-1)*eiko_m.z_spacing<<" m\n\n";
    
    if (reciprocity)
        std::cout<<"Reciprocity = True\n\n";
    else
        std::cout<<"Shots reciprocity = False\n\n";

    std::cout<<"Shot "<<shot_id+1<<" of "<<geometry[shots_type]->shots.all<<"\n";

    std::cout<<"Position (x,y,z) = ("<<geometry[shots_type]->shots.x[shot_id]<<", "
                                     <<geometry[shots_type]->shots.y[shot_id]<<", "
                                     <<geometry[shots_type]->shots.z[shot_id]<<") m\n\n";
}

void Eikonal::set_parameters(std::string file)
{
    eiko_m.x_samples = std::stoi(fm.catch_parameter("x_samples", file));
    eiko_m.y_samples = std::stoi(fm.catch_parameter("y_samples", file));
    eiko_m.z_samples = std::stoi(fm.catch_parameter("z_samples", file));

    eiko_m.total_samples = eiko_m.x_samples * eiko_m.y_samples * eiko_m.z_samples;

    eiko_m.x_spacing = std::stof(fm.catch_parameter("x_spacing", file));    
    eiko_m.y_spacing = std::stof(fm.catch_parameter("y_spacing", file));    
    eiko_m.z_spacing = std::stof(fm.catch_parameter("z_spacing", file));    

    slowness = new float[eiko_m.total_samples]();

    std::string vp_file = fm.catch_parameter("vp_file", file);
    fm.read_binary_float(vp_file, slowness, eiko_m.total_samples);
    
    for (int index = 0; index < eiko_m.total_samples; index++)
        slowness[index] = 1.0f / slowness[index];

    export_time_volume = fm.str2bool(fm.catch_parameter("export_time_volume", file));
    export_illumination = fm.str2bool(fm.catch_parameter("export_illumination",file));
    export_ray_position = fm.str2bool(fm.catch_parameter("export_ray_position", file));
    export_first_arrival = fm.str2bool(fm.catch_parameter("export_first_arrival", file));    

    ray_folder = fm.catch_parameter("ray_folder", file);
    time_volume_folder = fm.catch_parameter("time_volume_folder", file);
    illumination_folder = fm.catch_parameter("illumination_folder", file);
    first_arrival_folder = fm.catch_parameter("first_arrival_folder", file);

    shots_type = std::stoi(fm.catch_parameter("shot_geometry_type", file)); 
    nodes_type = std::stoi(fm.catch_parameter("node_geometry_type", file)); 
    reciprocity = fm.str2bool(fm.catch_parameter("reciprocity", file)); 

    geometry[0] = new Regular();
    geometry[1] = new Circular();

    geometry[shots_type]->set_parameters(file);    
    geometry[nodes_type]->set_parameters(file);

    geometry[shots_type]->export_positions(geometry[shots_type]->shots, "xyz_shots.txt");
    geometry[nodes_type]->export_positions(geometry[nodes_type]->nodes, "xyz_nodes.txt");
    
    if (reciprocity)
    {
        std::swap(geometry[shots_type]->shots.x, geometry[nodes_type]->nodes.x); 
        std::swap(geometry[shots_type]->shots.y, geometry[nodes_type]->nodes.y); 
        std::swap(geometry[shots_type]->shots.z, geometry[nodes_type]->nodes.z); 

        std::swap(geometry[shots_type]->shots.all, geometry[nodes_type]->nodes.all);
    }
    
    travel_time = new float[eiko_m.total_samples]();
    illumination = new float[eiko_m.total_samples]();
    first_arrival = new float[geometry[nodes_type]->nodes.all]();
}

void Eikonal::write_time_volume()
{
    if (export_time_volume)        
        fm.write_binary_float(time_volume_folder + "eikonal_nz" + std::to_string(eiko_m.z_samples) + "_nx" + std::to_string(eiko_m.x_samples) + "_ny" + std::to_string(eiko_m.y_samples) + "_shot_" + std::to_string(shot_id+1) + ".bin", travel_time, eiko_m.total_samples);
}

void Eikonal::write_illumination()
{
    if (export_illumination)
        fm.write_binary_float(illumination_folder + "illumination_nz" + std::to_string(eiko_m.z_samples) + "_nx" + std::to_string(eiko_m.x_samples) + "_ny" + std::to_string(eiko_m.y_samples) + ".bin", illumination, eiko_m.total_samples);
}

void::Eikonal::write_first_arrival()
{
    if (export_first_arrival) 
    {   
        int nx = eiko_m.x_samples;
        int nz = eiko_m.z_samples;

        float dx = eiko_m.x_spacing;
        float dy = eiko_m.y_spacing;
        float dz = eiko_m.z_spacing;

        for (int r = 0; r < geometry[nodes_type]->nodes.all; r++)
        {
            interpolate.x = geometry[nodes_type]->nodes.x[r];
            interpolate.y = geometry[nodes_type]->nodes.y[r];
            interpolate.z = geometry[nodes_type]->nodes.z[r];

            interpolate.x0 = floorf(interpolate.x / dx) * dx;
            interpolate.y0 = floorf(interpolate.y / dy) * dy;
            interpolate.z0 = floorf(interpolate.z / dz) * dz;

            interpolate.x1 = floorf(interpolate.x / dx) * dx + dx;
            interpolate.y1 = floorf(interpolate.y / dy) * dy + dy;
            interpolate.z1 = floorf(interpolate.z / dz) * dz + dz;

            int id = ((int)(interpolate.z / dz)) + 
                     ((int)(interpolate.x / dx))*nz + 
                     ((int)(interpolate.y / dy))*nx*nz;

            interpolate.c000 = travel_time[id];
            interpolate.c001 = travel_time[id + 1];
            interpolate.c100 = travel_time[id + nz]; 
            interpolate.c101 = travel_time[id + 1 + nz]; 
            interpolate.c010 = travel_time[id + nx*nz]; 
            interpolate.c011 = travel_time[id + 1 + nx*nz]; 
            interpolate.c110 = travel_time[id + nz + nx*nz]; 
            interpolate.c111 = travel_time[id + 1 + nz + nx*nz];

            first_arrival[r] = interpolate.trilinear();        
        }

        fm.write_binary_float(first_arrival_folder + "times_nr" + std::to_string(geometry[nodes_type]->nodes.all) + "_shot_" + std::to_string(shot_id+1) + ".bin", first_arrival, geometry[nodes_type]->nodes.all);
    }
}

void Eikonal::ray_tracing()
{
    if (export_ray_position || export_illumination)
    {        
        eiko_m.x_samples_b = eiko_m.x_samples + 2;    
        eiko_m.y_samples_b = eiko_m.y_samples + 2;    
        eiko_m.z_samples_b = eiko_m.z_samples + 2;    
    
        eiko_m.total_samples_b = eiko_m.x_samples_b * eiko_m.y_samples_b * eiko_m.z_samples_b;
    
        float * T = new float[eiko_m.total_samples_b];

        eiko_m.expand_fdm(travel_time, T);   

        int nxx = eiko_m.x_samples_b;
        int nzz = eiko_m.z_samples_b;

        float dx = eiko_m.x_spacing;
        float dy = eiko_m.y_spacing;
        float dz = eiko_m.z_spacing;
    
        float sx = geometry[shots_type]->shots.x[shot_id];
        float sy = geometry[shots_type]->shots.y[shot_id];
        float sz = geometry[shots_type]->shots.z[shot_id];

        int sIdx = (int)(sx / dx);
        int sIdy = (int)(sy / dy);
        int sIdz = (int)(sz / dz);

        int sId = sIdz + sIdx*eiko_m.z_samples + sIdy*eiko_m.x_samples*eiko_m.z_samples;     
        
        std::vector<float> xRay, yRay, zRay, iRay;

        int im, jm, km, rId;

        float rayStep = 0.2f * (dx + dy + dz) / 3.0f;

        for (int rayId = 0; rayId < geometry[nodes_type]->nodes.all; rayId++)
        {
            float zi = geometry[nodes_type]->nodes.z[rayId];
            float xi = geometry[nodes_type]->nodes.x[rayId];
            float yi = geometry[nodes_type]->nodes.y[rayId];

            im = (int)(zi / dz); 
            jm = (int)(xi / dx); 
            km = (int)(yi / dy); 

            rId = im + jm*eiko_m.z_samples + km*eiko_m.x_samples*eiko_m.z_samples;

            if (export_illumination) 
                illumination[rId] += rayStep;

            if (export_ray_position)
            {
                xRay.push_back(xi);
                yRay.push_back(yi);
                zRay.push_back(zi);
                iRay.push_back((float) rayId);
            }

            while (true)
            {
                int i = (int)(zi / dz) + 1;
                int j = (int)(xi / dx) + 1;
                int k = (int)(yi / dy) + 1;

                float dTz = (T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*dz);    
                float dTx = (T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*dx);    
                float dTy = (T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*dy);

                float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

                zi -= rayStep*dTz / norm; // z ray position update   
                xi -= rayStep*dTx / norm; // x ray position update   
                yi -= rayStep*dTy / norm; // y ray position update   

                im = (int)(zi / dz); 
                jm = (int)(xi / dx); 
                km = (int)(yi / dy); 

                rId = im + jm*eiko_m.z_samples + km*eiko_m.x_samples*eiko_m.z_samples;

                if (rId == sId)
                {
                    if (export_ray_position)
                    {
                        xRay.push_back(sx);
                        yRay.push_back(sy);
                        zRay.push_back(sz);

                        iRay.push_back((float) rayId);
                    }
                
                    if (export_illumination)
                    {
                        float finalDist = sqrtf(powf(sx - xi, 2.0f) + powf(sy - yi, 2.0f) + powf(sz - zi, 2.0f));        
                        illumination[rId] += finalDist;                
                    }            

                    break;
                }

                if (export_illumination) 
                    illumination[rId] += rayStep;

                if (export_ray_position)
                {
                    xRay.push_back(xi);
                    yRay.push_back(yi);
                    zRay.push_back(zi);
                    iRay.push_back((float) rayId);
                }
            }
        }

        if (export_ray_position)
        {
            std::ofstream fout;

            size_t allRayPoints = iRay.size();

            std::string xRayPath = ray_folder + "xRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shot_id+1) + ".bin";
            std::string yRayPath = ray_folder + "yRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shot_id+1) + ".bin";
            std::string zRayPath = ray_folder + "zRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shot_id+1) + ".bin";
            std::string iRayPath = ray_folder + "iRay_" + std::to_string(allRayPoints) + "_shot_" + std::to_string(shot_id+1) + ".bin";

            fout.open(xRayPath, std::ios::out);
            fout.write((char*)&xRay[0], xRay.size() * sizeof(xRay));
            std::cout<<"Binary file " + xRayPath<<" was succesfully written."<<std::endl;  
            fout.close();

            fout.open(yRayPath, std::ios::out);
            fout.write((char*)&yRay[0], yRay.size() * sizeof(yRay));
            std::cout<<"Binary file " + yRayPath<<" was succesfully written."<<std::endl;  
            fout.close();

            fout.open(zRayPath, std::ios::out);
            fout.write((char*)&zRay[0], zRay.size() * sizeof(zRay));
            std::cout<<"Binary file " + zRayPath<<" was succesfully written."<<std::endl;  
            fout.close();

            fout.open(iRayPath, std::ios::out);
            fout.write((char*)&iRay[0], iRay.size() * sizeof(iRay));
            std::cout<<"Binary file " + iRayPath<<" was succesfully written."<<std::endl;  
            fout.close();
        }

        std::vector<float>().swap(iRay);
        std::vector<float>().swap(xRay);
        std::vector<float>().swap(yRay);
        std::vector<float>().swap(zRay);

        delete[] T;
    }
}

