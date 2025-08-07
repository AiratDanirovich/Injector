#pragma once

#include <vector>
#include <iostream>

#include <Injector/Grids/Defines.h>

using VR = std::vector<GPN::RealType>;

void print_vector(const VR &data)
{
    for (auto i{0ull}; i < data.size(); ++i)
        std::cout << data[i] << "; ";
    std::cout << std::endl;
}