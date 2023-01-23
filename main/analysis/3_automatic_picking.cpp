# include <cmath>
# include <chrono>
# include <string>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "../../src/essentials/utils.hpp"

int main()
{
    auto utils = Utils();

    auto ti = std::chrono::system_clock::now();

    int nt = 3001;
    int traces = 100;
    int nodes_all = 441; 
    int shots_all = 10000;

    float dt = 1e-3f;
    float tlag = 0.10f;
    float tcut = 2.80f;

    int updated_nt = (int)(tcut/dt);

    float window = 0.080f;
    int iw = (int)(window/dt);

    int itime = (int)(tlag/dt);
    int ftime = (int)((tcut+tlag)/dt + 1);

    // loop of receiver gathers
    for (int node = 0; node < nodes_all; node++)
    {
        float * picks_all = new float[shots_all]();    
        float * seismic_all = utils.readBinaryFloat("../../inputs/seismograms/seismogram_" + std::to_string(nt) + "x" + std::to_string(shots_all) + "_shot_" + std::to_string(node+1) + ".bin", shots_all*nt);
        
        std::cout<<"Running node " + std::to_string(node+1) + " of " + std::to_string(nodes_all)<<std::endl;

        // loop of traces
        for (int trace = 0; trace < shots_all; trace++)
        {
            float * A = new float[updated_nt]();
            float * B = new float[updated_nt]();
            float * S = new float[updated_nt]();

            if (trace % (shots_all/10) == 0)
            {
                std::cout<<"    Running trace " + std::to_string(trace) + " of " + std::to_string(shots_all)<<std::endl;
            }

            // loop of time
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

            // Normalizing and picking the first break
            for (int k = 0; k < updated_nt; k++)
            {
                S[k] *= 1.0f / ampMax;
                
                if (S[k] > 1e-4f)
                {
                    picks_all[trace] = (k + iw + iw/2) * dt;            
                    break;
                }
            }    

            delete[] A;
            delete[] B;
            delete[] S;
        }

        utils.writeBinaryFloat("../../inputs/picks/rawPicks_shot_" + std::to_string(node+1) + "_" + std::to_string(shots_all) + "_samples.bin", picks_all, shots_all);

        delete[] picks_all;
        delete[] seismic_all;
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout << "\nRun time: " << elapsed_seconds.count() << " s\n";

    return 0;
}