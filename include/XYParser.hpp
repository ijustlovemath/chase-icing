#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
/* largely taken from https://stackoverflow.com/a/415553/2289030 */
template<typename U = double>
void load_xy_data(const std::string &filename, std::vector<U> &x, std::vector<U> &y)
{
    std::ifstream data(filename);

    std::string line;
    while(std::getline(data, line))
    {
        std::stringstream linestream(line);
        std::string cell;
        int i = 0;
        while(std::getline(linestream, cell, ','))
        {
            const U element = static_cast<U>(std::stod(cell));
            switch(i) {
            case 0:
                x.push_back(element);
                break;
            case 1:
                y.push_back(element);
                break;
            default: ;
            }
            ++i;
        }
        if(i < 2) { 
            // avoid mismatch if line is malformated
            x.pop_back();
            std::cerr << "problem loading data from " << filename << std::endl;
        }
    }
}
