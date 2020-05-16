#ifndef SERREFLECTION_TEST_STRUCTURES_1_H
#define SERREFLECTION_TEST_STRUCTURES_1_H

#include <serreflection/defines.hpp>

struct Inner_1
{
    int a;

    std::vector<double> b;

    std::string str;
};

SERIALIZIBLE_STRUCT(Inner_1, srfl::CheckModes::FORCE_SET,
                    (int, a, srfl::nan, 20, 30)(std::vector<double>, b, DEF_D())(std::string, str,
                                                                                 DEF_D()))

struct TestJSON_1
{
    int a;

    std::vector<double> b;

    double c;

    bool e;

    Inner_1 inner;

    std::vector<Inner_1> inners;
};

SERIALIZIBLE_STRUCT(TestJSON_1, srfl::CheckModes::FATAL,
                    (int, a, DEF_D())(std::vector<double>, b, srfl::nan, 4, 30)(double, c, DEF_D())(
                        bool, e, DEF_D())(Inner_1, inner, DEF_D())(std::vector<Inner_1>, inners,
                                                                   DEF_D()))

#endif
