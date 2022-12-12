# include <iostream>

# include "../../essentials/utils.hpp"
# include "../../essentials/model.hpp"
# include "../../essentials/geometry.hpp"
# include "../../eikonal/eikonal.hpp"
# include "../../tomography/tomography.hpp"

int main(int argc, char **argv)
{
    auto tomo = Tomography();

    tomo.setParameters(argv[1]);


    return 0;
}




