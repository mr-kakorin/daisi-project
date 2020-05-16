#ifndef SERREFLECTION_TEST_STRUCTURES_2_H
#define SERREFLECTION_TEST_STRUCTURES_2_H

#include <serreflection/defines.hpp>

namespace test_srfl
{
enum class En
{
    V1,
    V2
};

struct Inner_2
{
    int a;

    std::vector<double> b;

    std::string str;
};

struct TestJSON_2
{
    int a;

    std::vector<double> b;

    double c;

    bool e;

    Inner_2 inner;

    std::vector<Inner_2> inners;

    En v;
};
}

SERIALIZIBLE_ENUM(test_srfl::En, (V1)(V2))

SERIALIZIBLE_STRUCT(test_srfl::Inner_2, srfl::CheckModes::FORCE_SET,
                    (int, a, 20, 20, 30)(std::vector<double>, b, DEF_D())(std::string, str,
                                                                          DEF_D()))

SERIALIZIBLE_STRUCT(test_srfl::TestJSON_2, srfl::CheckModes::FORCE_SET,
                    (int, a, DEF_D())(std::vector<double>, b, srfl::nan, 4,
                                      30)(double, c, 1.4, -200, 200)(bool, e, 1, INF_D())(
                        test_srfl::Inner_2, inner, DEF_D())(std::vector<test_srfl::Inner_2>, inners,
                                                            DEF_D())(test_srfl::En, v, 1, INF_D()))

#endif
