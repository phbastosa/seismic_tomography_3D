# include <chrono>
# include <iostream>

# include "classic/classic.hpp"
# include "block_FIM/block_FIM.hpp"
# include "accurate_FSM/accurate_FSM.hpp"

# include "../utils/file_manager/file_manager.hpp"

int main(int argc, char **argv)
{
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    ti = std::chrono::system_clock::now();

    Eikonal * eikonal[] = 
    {
        new Classic(),
        new Block_FIM(),
        new Accurate_FSM()
    };

    int type = std::stoi(catch_parameter("eikonal_type", std::string(argv[1])));

    eikonal[type]->parameters = std::string(argv[1]);

    eikonal[type]->set_parameters();
    eikonal[type]->prepare_volumes();

    for (int i = 0; i < eikonal[type]->total_shots; i++)
    {
        eikonal[type]->shot_id = i;
        eikonal[type]->info_message();
        
        eikonal[type]->eikonal_equation();
        eikonal[type]->write_time_volume();
        eikonal[type]->write_first_arrival();
    }
    
    eikonal[type]->destroy_volumes();

    tf = std::chrono::system_clock::now();

    elapsed_seconds = tf - ti;
    std::cout<<"\nRun time = "<<elapsed_seconds.count()<<" s."<<std::endl;

    return 0;
}
