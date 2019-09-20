#include "Output.h"

void Output::openAll()
{
    for (size_t f = 0u; f < names.size(); ++f) {
        std::string filename = names[f] + ".dat";
        files.push_back(new std::ofstream());
        files.back()->open(filename, std::ios::binary);
        if (!files.back()->is_open()) {
            std::string msg = "Unable to open output file " + filename;
            throw std::runtime_error(msg);
        }
    }
}

void Output::closeAll()
{
    for (size_t f = 0u; f < names.size(); ++f) {
        files[f]->close();
    }
}
