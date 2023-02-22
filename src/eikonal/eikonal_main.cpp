# include <chrono>
# include <iostream>

# include "eikonal.hpp"
# include "classic/classic.hpp"
# include "block_FIM/block_FIM.hpp"
# include "accurate_FSM/accurate_FSM.hpp"

# include "../utils/file_manager/file_manager.hpp"

int main(int argc, char **argv)
{
    File_manager fm;

    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    ti = std::chrono::system_clock::now();

    Eikonal * eikonal[] = 
    {
        new Classic(),
        new Block_FIM(),
        new Accurate_FSM()
    };

    std::vector<std::string> formulation 
    {
        std::string("Podvin & Lecomte (1991)"),
        std::string("Jeong & Whitaker (2008)"),
        std::string("Noble, Gesret and Belayouni (2014)") 
    };

    int type = std::stoi(fm.catch_parameter("eikonal_type", std::string(argv[1])));

    eikonal[type]->set_parameters(std::string(argv[1]));
    
    eikonal[type]->prepare_volumes();

    for (eikonal[type]->shot_id = 0; eikonal[type]->shot_id < eikonal[type]->geometry[eikonal[type]->shots_type]->shots.all; eikonal[type]->shot_id++)
    {
        eikonal[type]->info_message();

        std::cout<<"Solving eikonal equation with the "<<formulation[type]<<" formulation\n";
        
        eikonal[type]->solve();
        eikonal[type]->write_time_volume();
        eikonal[type]->write_first_arrival();

        eikonal[type]->ray_tracing();
    }
    
    eikonal[type]->write_illumination();
    
    tf = std::chrono::system_clock::now();

    elapsed_seconds = tf - ti;
    std::cout<<"\nRun time = "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}
