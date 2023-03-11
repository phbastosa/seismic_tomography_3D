# include <cmath>
# include <chrono>
# include <string>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "../file_manager/file_manager.hpp"

int main(int argc, char**argv)
{
    auto ti = std::chrono::system_clock::now();

    int nt = 3001;
    int nodes_all = 11 * 16; 
    int shots_all = 140 * 100;

    float dt = 0.001f;
    
    float window = 0.01f;
    float amp_cut = 1e-15f;

    std::string data_folder = "../../../inputs/seismic_data/";    
    std::string pick_folder = "../../../inputs/picks/";

    int iw = (int)(window/dt)+1;

    for (int node = 0; node < nodes_all; node++)
    {
        float * picks_all = new float[shots_all]();    
        float * seismic_all = new float[shots_all * nt];
        
        read_binary_float(data_folder + "seismogram_" + std::to_string(nt) + "x" + std::to_string(shots_all) + "_shot_" + std::to_string(node+1) + ".bin", seismic_all, shots_all * nt);

        float seismicMaxAmp = *std::max_element(seismic_all, seismic_all + shots_all * nt);

        for (int i = 0; i < shots_all * nt; i++)
            seismic_all[i] *= 1.0f / seismicMaxAmp;

        std::cout<<"Running node " + std::to_string(node+1) + " of " + std::to_string(nodes_all)<<std::endl;

        for (int trace = 0; trace < shots_all; trace++)
        {
            float * A = new float[nt]();
            float * B = new float[nt]();
            float * S = new float[nt](); 

            if (trace % (shots_all/10) == 0)
            {
                std::cout<<"    Running trace " + std::to_string(trace) + " of " + std::to_string(shots_all)<<std::endl;
            }

            for (int time = iw; time < nt - iw; time++)
            {
                for (int k = 0; k < iw; k++)
                {
                    A[time] += seismic_all[(time - k) + trace*nt] + 1e-15f;
                    B[time] += seismic_all[(time + k) + trace*nt] + 1e-15f;
                }
                
                S[time] = fabsf(B[time]/A[time] * (B[time] - A[time]));
                S[time] *= (B[time] + A[time]) * (B[time] - A[time]) / powf(iw + time, 2.0f);
            }

            float ampMax = *std::max_element(S, S+nt);

            for (int k = 0; k < nt; k++)
            {
                S[k] *= 1.0f / ampMax;
                
                if (S[k] > amp_cut)
                {
                    picks_all[trace] = k * dt;            
                    break;
                }
            }    

            delete[] A;
            delete[] B;
            delete[] S;
        }

        write_binary_float(pick_folder + "rawPicks_shot_" + std::to_string(node+1) + "_" + std::to_string(shots_all) + "_samples.bin", picks_all, shots_all);

        delete[] picks_all;
        delete[] seismic_all;
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "\nRun time: " << elapsed_seconds.count() << " s\n";

    return 0;
}