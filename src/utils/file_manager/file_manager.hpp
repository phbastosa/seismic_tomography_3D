# ifndef FILE_MANAGER_CPP
# define FILE_MANAGER_CPP

# include <string>
# include <vector>

class File_manager
{
public:

    std::string parameter_file;

    bool str2bool(std::string s);

    void read_binary_float(std::string path, float * array, int n);
    void write_binary_float(std::string path, float *array, int n);

    void read_text_file(std::string path, std::vector<std::string> elements); 

    std::string catch_parameter(std::string target);

    std::vector<std::string> split(std::string s, char delimiter);
};

# endif