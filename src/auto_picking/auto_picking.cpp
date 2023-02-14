# include <cmath>
# include <chrono>
# include <string>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "../file_manager/file_manager.hpp"

int main()
{
    auto fm = File_manager();

    auto ti = std::chrono::system_clock::now();

    int nt = std::stoi(fm.catch_parameter("nt"));
    int nodes_all = std::stoi(fm.catch_parameter("nodes_all")); 
    int shots_all = std::stoi(fm.catch_parameter("shots_all"));

    float dt = std::stof(fm.catch_parameter("dt"));
    float tlag = std::stof(fm.catch_parameter("tlag"));
    float tcut = std::stof(fm.catch_parameter("tcut"));
    
    float window = std::stof(fm.catch_parameter("pick_window"));
    float amp_cut = std::stof(fm.catch_parameter("amp_cut"));
    float pick_lag = std::stof(fm.catch_parameter("pick_lag"));

    std::string data_folder = fm.catch_parameter("data_folder");    
    std::string pick_folder = fm.catch_parameter("pick_folder");

    int updated_nt = (int)(tcut/dt)+1;

    int iw = (int)(window/dt)+1;

    int itime = (int)(tlag/dt)+1;
    int ftime = (int)((tcut+tlag)/dt)+1;

    for (int node = 0; node < nodes_all; node++)
    {
        float * picks_all = new float[shots_all]();    
        float * seismic_all = new float[shots_all * nt];
        
        fm.read_binary_float(data_folder + "seismogram_" + std::to_string(nt) + "x" + std::to_string(shots_all) + "_shot_" + std::to_string(node+1) + ".bin", seismic_all, shots_all * nt);
        
        std::cout<<"Running node " + std::to_string(node+1) + " of " + std::to_string(nodes_all)<<std::endl;

        for (int trace = 0; trace < shots_all; trace++)
        {
            float * A = new float[updated_nt]();
            float * B = new float[updated_nt]();
            float * S = new float[updated_nt]();

            if (trace % (shots_all/10) == 0)
            {
                std::cout<<"    Running trace " + std::to_string(trace) + " of " + std::to_string(shots_all)<<std::endl;
            }

            for (int time = iw; time < updated_nt - iw; time++)
            {
                for (int k = 0; k < iw; k++)
                {
                    A[time] += seismic_all[(itime+time-k) + trace*nt] + 1e-15f;
                    B[time] += seismic_all[(itime+time+k) + trace*nt] + 1e-15f;
                }
                
                S[time] = fabsf(B[time]/A[time] * (B[time] - A[time]));
                S[time] *= (B[time] + A[time]) * (B[time] - A[time]) / powf(iw + time, 2.0f);
            }

            float ampMax = *std::max_element(S, S+updated_nt);

            for (int k = 0; k < updated_nt; k++)
            {
                S[k] *= 1.0f / ampMax;
                
                if (S[k] > amp_cut)
                {
                    picks_all[trace] = (k + pick_lag*iw) * dt;            
                    break;
                }
            }    

            delete[] A;
            delete[] B;
            delete[] S;
        }

        fm.write_binary_float(pick_folder + "rawPicks_shot_" + std::to_string(node+1) + "_" + std::to_string(shots_all) + "_samples.bin", picks_all, shots_all);

        delete[] picks_all;
        delete[] seismic_all;
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "\nRun time: " << elapsed_seconds.count() << " s\n";

    return 0;
}