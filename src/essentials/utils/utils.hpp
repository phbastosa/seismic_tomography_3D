# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>
# include <string>

class Utils
{
private:

public:
    
    typedef struct 
    {
        float x;
        float y; 
        
    } point2D;

    static bool str2bool(std::string s);

    static std::vector<float> linspace(float start, float end, int num);        

    static std::vector<std::string> split(std::string s, char delimiter);
};

# endif
