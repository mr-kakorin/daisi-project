#ifndef SERREFLECTION_WRITE_JSON_H
#define SERREFLECTION_WRITE_JSON_H

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <boostlog0/boostlog.h>

#include <serreflection/private/write_json_private.hpp>

namespace srfl
{
template <typename T>
std::unique_ptr<boost::property_tree::ptree> write_json_pt(const T& object) noexcept
{
    try
    {
        return write_json_p(object);
    }
    catch (const std::exception& ex)
    {
        LOG(sev_lvl::error) << ex.what();
        return nullptr;
    }
}

template <class T>
std::string write_json_str(const T& object) noexcept
{
    auto pt = write_json_pt(object);

    if (!pt)
    {
        return "";
    }

    std::ostringstream oss;
    boost::property_tree::json_parser::write_json(oss, *pt);

    return oss.str();
}
}

#endif
