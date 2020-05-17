// copy of propagate_const adopted for gcc 4.9

#ifndef COMM_TOOLS_PROPAGATE_CONST_HPP
#define COMM_TOOLS_PROPAGATE_CONST_HPP

#ifdef __GNUC__
#if __GNUC__ < 6 
#include <common_tools/propagate_const_.hpp>
#else
#include <experimental/propagate_const>
#endif

namespace commtools
{

template <class T>
using pc_unique_ptr = std::experimental::propagate_const<std::unique_ptr<T>>;
template <class T>
using pc_shared_ptr = std::experimental::propagate_const<std::shared_ptr<T>>;

}

#else
#include <memory>
namespace commtools
{
	template <class T>
	using pc_unique_ptr = std::unique_ptr<T>;
	template <class T>
	using pc_shared_ptr = std::shared_ptr<T>;
}
#endif

#endif
