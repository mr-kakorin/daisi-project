
#ifndef COMMON_TOOLS_JSON_HELPER_H
#define COMMON_TOOLS_JSON_HELPER_H

#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace commtools
{
template <typename T>
std::vector<T> json_vector(const boost::property_tree::ptree& pt, const std::string& key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
    {
        r.push_back(item.second.get_value<T>());
    }
    return r;
}
bool read_property_tree(const boost::property_tree::ptree& config,
                        const std::string&                 error_msg_prefix,
                        const std::function<bool(const boost::property_tree::ptree&)>& operation);

bool read_property_tree_string(
    const std::string& config, const std::string& error_msg_prefix,
    const std::function<bool(const boost::property_tree::ptree&)>& operation);

bool read_property_tree_file(
    const std::string& path, const std::string& error_msg_prefix,
    const std::function<bool(const boost::property_tree::ptree&)>& operation);

std::shared_ptr<boost::property_tree::ptree> readJSONString(const std::string& config);

std::shared_ptr<boost::property_tree::ptree> readJSONFile(const std::string& path);
}

#endif // COMMON_TOOLS_JSON_HELPER_H