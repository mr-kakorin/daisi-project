
#ifndef SERREFLECTION_WRITE_JSON_PRIVATE_H
#define SERREFLECTION_WRITE_JSON_PRIVATE_H

#include <serreflection/private/json_writer.hpp>

#include <serreflection/private/universal_reader_writer.hpp>

namespace srfl
{
template <typename T>
std::unique_ptr<boost::property_tree::ptree> write_json_p(const T& object)
{
    return read_write<JSONWriter<T>, T, boost::property_tree::ptree>(object);
}
}

#endif
