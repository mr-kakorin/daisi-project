#include <algorithm>
#include <cstdlib>

#include "tools.h"

namespace notk
{
namespace tools
{
std::vector<size_t> get_rand_indexes(const size_t length)
{
    std::vector<size_t> indexes(length);
    std::vector<size_t> result(length);

    for (size_t dim = 0; dim < length; dim++)
    {
        indexes[dim] = dim;
    }

    for (size_t dim = 0; dim < length; dim++)
    {
        size_t ind   = rand() % (indexes.size());
        size_t value = indexes[ind];
        result[dim]  = value;
        indexes.erase(indexes.begin() + ind);
    }
    return result;
}
std::vector<std::vector<size_t>> get_array_of_rand_indexes(const size_t length, const size_t n)
{
    size_t real_size = n;
    if (length < 20)
    {
        real_size = std::min(n, factorial(length));
    }

    std::vector<std::vector<size_t>> result(real_size);

    result[0] = get_rand_indexes(length);
    for (size_t k = 1; k < real_size; k++)
    {
        while (true)
        {
            result[k] = get_rand_indexes(length);

            bool flagBreak = true;

            for (size_t kk = 0; kk < k; kk++)
            {
                if (std::equal(result[kk].begin(), result[kk].end(), result[k].begin()))
                {
                    flagBreak = false;
                    break;
                }
            }
            if (flagBreak)
            {
                break;
            }
        }
    }
    return result;
}
}
}