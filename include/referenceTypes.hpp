#ifndef REFERENCETYPES_HPP
#define REFERENCETYPES_HPP 

#include <functional>
#include <array>
#include <map>


template<typename T>
using IJ_map =  std::map<std::array<int,2>,T>; //A square map

template<typename T>
using I_map= std::map<int, T>; //equivalent to a vector with (possible) missing entries

#endif