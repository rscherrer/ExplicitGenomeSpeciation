#include "doMain.h"
#include <iostream>
#include <vector>
#include <string>


/// Program to run the main function
int doMain(const std::vector<std::string> &args)
{

    try
    {
        // Should return an error if there are more than one argument
        if (args.size() > 2u)
        {
            throw std::runtime_error("More than one argument was supplied");
        }
    }
    catch (const std::runtime_error &err)
    {
        std::cout << "Exception:" << err.what() << '\n';
        return 1;
    }

    return 0;
}
