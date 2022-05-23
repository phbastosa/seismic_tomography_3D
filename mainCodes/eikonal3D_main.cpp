# include "../src/essentials/inout/inout.hpp"
# include "../src/essentials/model/model.hpp"
# include "../src/essentials/utils/utils.hpp"
# include "../src/essentials/geometry/geometry.hpp"

# include "../src/simulations/eikonal/eikonal.hpp"

int main(int argc, char **argv)
{
    auto eikonal = Eikonal3D(argv);

    for (int i = 0; i < eikonal.g3D.ns; i++)
    {
        eikonal.shotId = i;

        if (eikonal.eikonalType)
        {
            eikonal.podvin3D();
        }
        else 
        {
            eikonal.fim3D();
        }
    }    

    eikonal.deleteVolumes();

    return 0;
}
