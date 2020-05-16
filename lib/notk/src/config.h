#ifndef NOTK_CONFIG_H
#define NOTK_CONFIG_H

#include <serreflection/defines.hpp>

#include "enums.h"

namespace notk
{
struct NOTKConfig
{
    bool allow_maximization;
    int  log_each_iteration;

    std::string version;

    InitialBordersType borders_type;
};
struct NOTKBorders
{
    int problem_dimension;

    std::vector<double> x_left;
    std::vector<double> x_right;
    std::vector<double> x_0;
};
}

SERIALIZIBLE_STRUCT(notk::NOTKConfig, srfl::CheckModes::FATAL,
                    (bool, allow_maximization, DEF_D())(int, log_each_iteration, srfl::nan, 1,
                                                        srfl::inf)(std::string, version, DEF_D())(
                        std::string, version, DEF_D())(notk::InitialBordersType, borders_type,
                                                       DEF_D()))

SERIALIZIBLE_STRUCT(notk::NOTKBorders, srfl::CheckModes::FATAL,
                    (int, problem_dimension, 1, 1, srfl::inf)(std::vector<double>, x_left, DEF_D())(
                        std::vector<double>, x_right, DEF_D())(std::vector<double>, x_0, DEF_D()))

#endif
