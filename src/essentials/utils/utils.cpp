# include <vector>
# include <string>
# include <sstream>
# include <algorithm>

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

std::vector<std::string> Utils::split(std::string s, char delimiter)
{
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) 
        tokens.push_back(token);
   
    return tokens;
}

bool Utils::str2bool(std::string s)
{
    bool b;

    std::for_each(s.begin(), s.end(), [](char & c) {c = ::tolower(c);});
    std::istringstream(s) >> std::boolalpha >> b;

    return b;
}
