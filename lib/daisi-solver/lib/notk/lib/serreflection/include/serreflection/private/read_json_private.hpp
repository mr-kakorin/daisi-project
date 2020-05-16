#ifndef SERREFLECTION_READ_JSON_PRIVATE_H
#define SERREFLECTION_READ_JSON_PRIVATE_H

#include <serreflection/private/json_reader.hpp>

#include <serreflection/private/universal_reader_writer.hpp>

namespace srfl
{
template <typename T, typename Tpt>
std::unique_ptr<T> read_json_p(Tpt&& pt)
{
    return read_write<JSONReader<T>, T, T>(std::forward<Tpt>(pt));
}
template <typename T, typename Tpt>
void read_json_p(Tpt&& pt, T& object)
{
    JSONReader<T> reader(std::forward<Tpt>(pt), object);

    using range = boost::mpl::range_c<size_t, 0, boost::fusion::result_of::size<T>::value>;

    boost::fusion::for_each(range(), reader);

    reader.get().release();
}
}

#endif
