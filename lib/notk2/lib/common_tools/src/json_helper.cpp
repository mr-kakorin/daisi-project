
#include <boostlog0/boostlog.h>

#include <common_tools/json_helper.h>

namespace commtools
{
bool read_property_tree(const boost::property_tree::ptree& config,
                        const std::string&                 error_msg_prefix,
                        const std::function<bool(const boost::property_tree::ptree&)>& operation)
{
    try
    {
        return operation(config);
    }

    catch (const boost::property_tree::json_parser_error& e)
    {
        LOG(sev_lvl::error) << error_msg_prefix << e.message();
    }
    catch (const std::exception& e)
    {
        LOG(sev_lvl::error) << error_msg_prefix << e.what();
    }
    catch (...)
    {
        LOG(sev_lvl::error) << error_msg_prefix << ": unknown error";
    }
    return false;
}

std::shared_ptr<boost::property_tree::ptree> readJSONString(const std::string& config)
{
    std::stringstream ss;
    ss << config;

    auto pt = std::make_shared<boost::property_tree::ptree>();
    try
    {
        boost::property_tree::read_json(ss, *pt);
        return pt;
    }
    catch (const boost::property_tree::json_parser_error& e)
    {
        LOG(sev_lvl::error) << "readJSONString::error: " << e.message();
    }
    catch (const std::exception& e)
    {
        LOG(sev_lvl::error) << "readJSONFile::error: " << e.what();
    }
    catch (...)
    {
        LOG(sev_lvl::error) << "readJSONString::unknown error";
    }
    return nullptr;
}
std::shared_ptr<boost::property_tree::ptree> readJSONFile(const std::string& path)
{

    auto pt = std::make_shared<boost::property_tree::ptree>();
    try
    {
        boost::property_tree::read_json(path, *pt);
        return pt;
    }
    catch (const boost::property_tree::json_parser_error& e)
    {
        LOG(sev_lvl::error) << "readJSONFile::error: " << e.message();
    }
    catch (const std::exception& e)
    {
        LOG(sev_lvl::error) << "readJSONFile::error: " << e.what();
    }
    catch (...)
    {
        LOG(sev_lvl::error) << "readJSONFile::unknown error";
    }
    return nullptr;
}
bool read_property_tree_string(
    const std::string& config, const std::string& error_msg_prefix,
    const std::function<bool(const boost::property_tree::ptree&)>& operation)
{
    auto pt = readJSONString(config);
    return !pt ? false : read_property_tree(*pt, error_msg_prefix, operation);
}

bool read_property_tree_file(
    const std::string& path, const std::string& error_msg_prefix,
    const std::function<bool(const boost::property_tree::ptree&)>& operation)
{
    auto pt = readJSONFile(path);
    return !pt ? false : read_property_tree(*pt, error_msg_prefix, operation);
}
}