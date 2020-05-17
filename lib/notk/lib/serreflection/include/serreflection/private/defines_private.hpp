#ifndef SERREFLECTION_DEFINES_PRIVATE_H
#define SERREFLECTION_DEFINES_PRIVATE_H

#include <array>
#include <map>

#include <boost/fusion/adapted/struct/adapt_struct.hpp>

#include <boostlog0/boostlog.h>

#define FILLER_X4(x1, x2, x3, x4, x5) ((x1, x2, x3, x4, x5)) FILLER_Y4
#define FILLER_Y4(x1, x2, x3, x4, x5) ((x1, x2, x3, x4, x5)) FILLER_X4
#define FILLER_X40
#define FILLER_Y40

#define FILLER_X1(x1) ((x1)) FILLER_Y1
#define FILLER_Y1(x1) ((x1)) FILLER_X1
#define FILLER_X10
#define FILLER_Y10

#define TO_MAP_INNER(R, vals, i, tuple)                                          \
    BOOST_PP_COMMA_IF(i)                                                         \
    {                                                                            \
        BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(5, 1, tuple))                     \
        , std::array<double, 3>                                                  \
        {                                                                        \
            BOOST_PP_TUPLE_ELEM(5, 2, tuple)                                     \
            , BOOST_PP_TUPLE_ELEM(5, 3, tuple), BOOST_PP_TUPLE_ELEM(5, 4, tuple) \
        }                                                                        \
    }

#define TO_MAP(seq)                                   \
    {                                                 \
        BOOST_PP_SEQ_FOR_EACH_I(TO_MAP_INNER, ~, seq) \
    }

#define TO_BOOST_FUSION(R, vals, i, tuple) \
    (BOOST_PP_TUPLE_ELEM(5, 0, tuple), BOOST_PP_TUPLE_ELEM(5, 1, tuple))

#define SERIALIZIBLE_STRUCT_BASE(NAME, MODE, ATTRIBUTES)                                           \
    BOOST_FUSION_ADAPT_STRUCT(NAME, BOOST_PP_SEQ_FOR_EACH_I(TO_BOOST_FUSION, ~, ATTRIBUTES))       \
    namespace srfl                                                                                 \
    {                                                                                              \
        namespace                                                                                  \
        {                                                                                          \
            template <>                                                                            \
            struct StructChecker<NAME>                                                             \
            {                                                                                      \
                const static std::map<std::string, std::array<double, 3>> ranges;                  \
                const static srfl::CheckModes mode = MODE;                                         \
                const static std::string      class_name;                                          \
                template <typename Tval>                                                           \
                static bool check_value(                                                           \
                    Tval &value, const std::string &name,                                          \
                    typename std::enable_if<std::is_fundamental<Tval>::value>::type * = 0)         \
                {                                                                                  \
                    if (srfl::CheckModes::NONE == mode)                                            \
                    {                                                                              \
                        return true;                                                               \
                    }                                                                              \
                    auto it = ranges.find(name);                                                   \
                    if (ranges.end() == it)                                                        \
                    {                                                                              \
                        throw std::runtime_error("check_value::Error : unexpected variable name"); \
                    }                                                                              \
                    auto left = it->second[1];                                                     \
                    auto right = it->second[2];                                                    \
                    if (!(value >= left && value <= right))                                        \
                    {                                                                              \
                        auto force_val = value < left ? left : right;                              \
                        auto msg = "Variable: " + class_name + "::" + name + " = " +               \
                                   std::to_string(value) + " is out of expected range [" +         \
                                   std::to_string(left) + ", " + std::to_string(right) + "] ";     \
                        switch (mode)                                                              \
                        {                                                                          \
                            case srfl::CheckModes::FATAL:                                          \
                                throw std::runtime_error(msg);                                     \
                            case srfl::CheckModes::FORCE_SET:                                      \
                                LOG(sev_lvl::warning)                                              \
                                    << msg << "Force set its value to " << force_val;              \
                                value = force_val;                                                 \
                                return true;                                                       \
                            case srfl::CheckModes::WARNING:                                        \
                                LOG(sev_lvl::warning) << msg;                                      \
                                return true;                                                       \
                            default:                                                               \
                                throw std::runtime_error("check_value::Error : unexpected mode");  \
                        }                                                                          \
                    }                                                                              \
                    return true;                                                                   \
                }                                                                                  \
                template <typename Tvector>                                                        \
                static bool check_value(                                                           \
                    std::vector<Tvector> &element, const std::string &name,                        \
                    typename std::enable_if<std::is_fundamental<Tvector>::value>::type * = 0)      \
                {                                                                                  \
                    bool result = true;                                                            \
                    std::for_each(std::begin(element), std::end(element),                          \
                                  [&result, &name](Tvector &value) {                               \
                                      if (!check_value(value, name))                               \
                                      {                                                            \
                                          result = false;                                          \
                                      }                                                            \
                                  });                                                              \
                    return result;                                                                 \
                }                                                                                  \
                static double get_default_value(const std::string &name)                           \
                {                                                                                  \
                    auto it = ranges.find(name);                                                   \
                    if (ranges.end() == it)                                                        \
                    {                                                                              \
                        throw std::runtime_error(                                                  \
                            "get_default_value::Error : unexpected variable name");                \
                    }                                                                              \
                    return it->second[0];                                                          \
                }                                                                                  \
            };                                                                                     \
            const std::map<std::string, std::array<double, 3>> StructChecker<NAME>::ranges =       \
                TO_MAP(ATTRIBUTES);                                                                \
            const std::string StructChecker<NAME>::class_name = #NAME;                             \
        }                                                                                          \
    }

#define EXPAND_SEQ(R, vals, i, tuple) \
    BOOST_PP_COMMA_IF(i)              \
    BOOST_PP_TUPLE_ELEM(1, 0, tuple)

#define ENUM_TO_MAP_INNER(R, NAME, i, tuple)                 \
    BOOST_PP_COMMA_IF(i)                                     \
    {                                                        \
        BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(1, 0, tuple)) \
        , NAME::BOOST_PP_TUPLE_ELEM(1, 0, tuple)             \
    }

#define ENUM_TO_MAP_INNER_INV(R, NAME, i, tuple)                 \
    BOOST_PP_COMMA_IF(i)                                         \
    {                                                            \
        NAME::BOOST_PP_TUPLE_ELEM(1, 0, tuple),                  \
            BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(1, 0, tuple)) \
    }

#define ENUM_TO_MAP(seq, NAME)                                \
    {                                                         \
        BOOST_PP_SEQ_FOR_EACH_I(ENUM_TO_MAP_INNER, NAME, seq) \
    }

#define ENUM_TO_MAP_INV(seq, NAME)                                \
    {                                                             \
        BOOST_PP_SEQ_FOR_EACH_I(ENUM_TO_MAP_INNER_INV, NAME, seq) \
    }

#define SERIALIZIBLE_ENUM_MAP_BASE(NAME, ATTRIBUTES)                            \
    namespace srfl                                                              \
    {                                                                           \
        namespace                                                               \
        {                                                                       \
            template <>                                                         \
            struct EnumParser<NAME>                                             \
            {                                                                   \
                const static std::map<std::string, NAME> parse_map;             \
                const static std::map<NAME, std::string> inv_parse_map;         \
            };                                                                  \
            const std::map<std::string, NAME> EnumParser<NAME>::parse_map =     \
                ENUM_TO_MAP(ATTRIBUTES, NAME);                                  \
            const std::map<NAME, std::string> EnumParser<NAME>::inv_parse_map = \
                ENUM_TO_MAP_INV(ATTRIBUTES, NAME);                              \
        }                                                                       \
    }

#define SERIALIZIBLE_ENUM_BUILD_BASE(NAME, ATTRIBUTES)     \
    enum class NAME                                        \
    {                                                      \
        BOOST_PP_SEQ_FOR_EACH_I(EXPAND_SEQ, ~, ATTRIBUTES) \
    };                                                     \
    SERIALIZIBLE_ENUM_MAP_BASE(NAME, ATTRIBUTES)

namespace srfl
{
    namespace
    {
        template <typename T>
        struct EnumParser
        {
        };
    }
    namespace
    {
        template <typename T>
        struct StructChecker
        {
        };
    }

    template <typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args &&... args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
}

#endif
