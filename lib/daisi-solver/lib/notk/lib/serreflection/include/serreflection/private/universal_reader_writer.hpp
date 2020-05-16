#ifndef SERREFLECTION_UNIVERSAL_READER_WRITER_H
#define SERREFLECTION_UNIVERSAL_READER_WRITER_H

#include <memory>

namespace srfl
{
template <class TMan, class TObj, typename TRes, typename TSrc>
std::unique_ptr<TRes> read_write(TSrc&& src)
{
    TMan manipulator(std::forward<TSrc>(src));

    using range = boost::mpl::range_c<size_t, 0, boost::fusion::result_of::size<TObj>::value>;

    boost::fusion::for_each(range(), manipulator);

    return std::move(manipulator.get());
}
}

#endif
