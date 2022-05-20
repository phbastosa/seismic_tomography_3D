# include "utils.hpp"

# include "../inout/inout.hpp"

int main(int argc, char **argv)
{
    std::vector<float> x = Utils::linspace(0.25f,0.75f,5);

    for (int i = 0; i < x.size(); i++) std::cout<<x[i]<<std::endl;

    Utils::point2D p;

    p.x = 20.4f;
    p.y = 50.3f;

    std::cout<<"Ponto ("<<p.x<<", "<<p.y<<") correctly written!"<<std::endl;

    return 0;
}
