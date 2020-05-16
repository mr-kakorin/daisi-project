#ifndef NOTK_OPT_TYPES_H
#define NOTK_OPT_TYPES_H

#include <utility>
#include <vector>

namespace notk
{
template <class Targ>
using borders_t = std::pair<std::vector<Targ>, std::vector<Targ>>;

template <class Targ, class Tfit>
using it_res_t = std::pair<std::vector<Targ>, Tfit>;
}

#endif
