# include <vector>

# include "utils.hpp"

std::vector<float> Utils::linspace(float start, float end, int num)
{
    std::vector<float> linspaced;
    
    if (num == 0) return linspaced;
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    } 

    linspaced.reserve(num);

    float delta = (end - start) / (num - 1);

    for (int i = 0; i < num; i++)
    {
        linspaced.emplace_back(start + (float)(delta*i));
    }

    return linspaced;
}
