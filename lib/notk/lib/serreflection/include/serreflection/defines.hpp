#ifndef SERREFLECTION_DEFINES_H
#define SERREFLECTION_DEFINES_H

#include <limits>

#include <serreflection/private/defines_private.hpp>

namespace srfl
{
    enum class CheckModes
    {
        FATAL,      // return nullptr
        FORCE_SET,  // force set value to nearest border
        WARNING,    // just print warning message
        NONE        // ignore
    };

    const static double inf = std::numeric_limits<double>::infinity();
    const static double nan = std::numeric_limits<double>::quiet_NaN();
}

#define SERIALIZIBLE_STRUCT(NAME, MODE, ATTRIBUTES) \
    SERIALIZIBLE_STRUCT_BASE(NAME, MODE, BOOST_PP_CAT(FILLER_X4 ATTRIBUTES, 0))

#define SERIALIZIBLE_ENUM_BUILD(NAME, ATTRIBUTES) \
    SERIALIZIBLE_ENUM_BUILD_BASE(NAME, BOOST_PP_CAT(FILLER_X1 ATTRIBUTES, 0))

#define SERIALIZIBLE_ENUM(NAME, ATTRIBUTES) \
    SERIALIZIBLE_ENUM_MAP_BASE(NAME, BOOST_PP_CAT(FILLER_X1 ATTRIBUTES, 0))

#define INF_D() -srfl::inf, srfl::inf
#define DEF_D() srfl::nan, INF_D()

#endif
