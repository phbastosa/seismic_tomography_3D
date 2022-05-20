# ifndef UTILS_HPP
# define UTILS_HPP

# include <vector>

class Utils
{
private:

public:
    
    typedef struct 
    {
        float x;
        float y; 
        
    } point2D;

    static std::vector<float> linspace(float start, float end, int num);        
};

# endif
