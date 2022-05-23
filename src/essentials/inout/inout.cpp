# include <string>
# include <fstream>
# include <iostream>
# include <algorithm>

# include "inout.hpp"

void InOut::readBinaryFloat(std::string path, float *array, int n)
{
    std::ifstream file(path, std::ios::binary);
    file.read((char *) array, n * sizeof(float));
    file.close();    
}

void InOut::writeBinaryFloat(std::string path, float *array, int n)
{
    std::ofstream file(path, std::ios::binary);
    file.write((char *) array, n * sizeof(float));
    file.close();
}

std::string InOut::catchParameter(std::string target, std::string file)
{
    char spaces = ' ';
    char comment = '#';

    std::string line;
    std::string variable;

    std::ifstream parameters(file);

    if (parameters.is_open())
    {
        while (getline(parameters, line))
        {           
            if ((line.front() != comment) && (line.front() != spaces))        
            {
                if (line.find(target) == 0)
                {
                    for (int i = line.find("=")+2; i < line.size(); i++)
                    {    
                        if (line[i] == '#') break;
                        variable += line[i];            
                    }

                    break;
                }
            }                 
        }
        parameters.close();
    }        
    else std::cout<<"Unable to open a parameter file!"<<std::endl;

    // Quality control for file paths

    if (variable.find('"') == 0)
    {
        remove(variable.begin(), variable.end(), '"');
    }
    else if (variable.find("[") == 0)
    {
        remove(variable.begin(), variable.end(), '[');
        remove(variable.begin(), variable.end(), ']');
    }

    variable.erase(remove(variable.begin(), variable.end(), ' '), variable.end());

    return variable;
}
