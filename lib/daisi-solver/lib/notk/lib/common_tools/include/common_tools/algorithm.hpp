
#ifndef COMMON_TOOLS_ALGORITHM_H
#define COMMON_TOOLS_ALGORITHM_H

#include <boost/algorithm/string.hpp>

namespace commtools
{
template <class T>
T strsplit_adapter(const std::string& data, const std::string& split)
{
    T result;

    boost::algorithm::split(result, data, boost::is_any_of(split));
    result.erase(
        std::remove_if(result.begin(), result.end(), [](const auto& val) { return val.empty(); }),
        result.end());

    return result;
}
}

#endif
