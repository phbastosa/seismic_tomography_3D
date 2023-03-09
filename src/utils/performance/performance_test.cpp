# include <chrono>
# include <fstream>
# include <iostream>

# include "../../eikonal/classic/classic.hpp"
# include "../../eikonal/block_FIM/block_FIM.hpp"
# include "../../eikonal/accurate_FSM/accurate_FSM.hpp"

# include "../file_manager/file_manager.hpp"

int main(int argc, char **argv)
{
    std::chrono::duration<double> elapsed_seconds;
    std::chrono::_V2::system_clock::time_point ti, tf;

    Eikonal * eikonal[] = 
    {
        new Classic(),
        new Block_FIM(),
        new Accurate_FSM()
    };

    std::ofstream file("eikonal_run_time.txt", std::ios::app);

    std::vector<std::string> method {"Podvin & Lecomte (1991)           ", 
                                     "Jeong & Whitaker (2008)           ", 
                                     "Noble, Gesret and Belayouni (2014)"};

    for (int type = 0; type < 3; type++) 
    {
        eikonal[type]->parameters = std::string(argv[1]);

        eikonal[type]->set_parameters();
        eikonal[type]->prepare_volumes();

        ti = std::chrono::system_clock::now();

        for (int i = 0; i < eikonal[type]->total_shots; i++)
        {
            eikonal[type]->shot_id = i;
            eikonal[type]->info_message();
            
            eikonal[type]->eikonal_equation();
            eikonal[type]->write_time_volume();
            eikonal[type]->write_first_arrival();
        }

        tf = std::chrono::system_clock::now();

        elapsed_seconds = tf - ti;

        file <<method[type]<<" = "<< elapsed_seconds.count() << "s.\n";

        eikonal[type]->destroy_volumes();
    }

    return 0;
}
