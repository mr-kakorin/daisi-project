#ifndef SERREFLECTION_READ_JSON_H
#define SERREFLECTION_READ_JSON_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <boostlog0/boostlog.h>

#include <serreflection/private/read_json_private.hpp>

namespace srfl
{
template <typename T, typename Tpt>
bool read_json(Tpt&& pt, T& object) noexcept
{
    try
    {
        read_json_p(std::forward<Tpt>(pt), object);
        return true;
    }
    catch (const std::exception& ex)
    {
        LOG(sev_lvl::error) << ex.what();
        return false;
    }
}

template <typename T, typename Tpt>
std::unique_ptr<T> read_json(Tpt&& pt) noexcept
{
    try
    {
        return read_json_p<T>(std::forward<Tpt>(pt));
    }
    catch (const std::exception& ex)
    {
        LOG(sev_lvl::error) << ex.what();
        return nullptr;
    }
}

template <class T>
std::unique_ptr<T> read_json_string(const std::string& json) noexcept
{
    std::stringstream ss;
    ss << json;
    boost::property_tree::ptree pt;

    try
    {
        boost::property_tree::read_json(ss, pt);
    }
    catch (const std::exception& ex)
    {
        LOG(sev_lvl::error) << ex.what();
        return nullptr;
    }

    return read_json<T, boost::property_tree::ptree&&>(std::move(pt));
}
}

#endif
