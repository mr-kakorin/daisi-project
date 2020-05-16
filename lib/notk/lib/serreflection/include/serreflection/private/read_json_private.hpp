#ifndef SERREFLECTION_READ_JSON_PRIVATE_H
#define SERREFLECTION_READ_JSON_PRIVATE_H

#include <serreflection/private/json_reader.hpp>

#include <serreflection/private/universal_reader_writer.hpp>

namespace srfl
{
    template <typename T, typename Tpt>
    std::unique_ptr<T> read_json_p(Tpt &&pt)
    {
        return read_write<JSONReader<T>, T, T>(std::forward<Tpt>(pt));
    }
}

#endif
