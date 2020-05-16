#ifndef NOTK_TOOLS_H
#define NOTK_TOOLS_H

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

namespace notk
{
namespace tools
{

template <class Tfit, class Targ>
std::vector<Tfit> calc_grad(std::function<Tfit(const std::vector<Targ>&)> function,
                            const std::vector<Targ>& arg, const double rel_step)
{
    std::vector<Tfit> result(arg.size());

    for (size_t i = 0; i < result.size(); i++)
    {
        auto tmp_arg = arg;

        double dx;

        if (std::abs(tmp_arg[i]) > 1e-15)
        {
            dx = std::abs(tmp_arg[i] * rel_step);
        }
        else
        {
            dx = rel_step * rel_step;
        }
        tmp_arg[i] = tmp_arg[i] + dx;

        result[i] = (function(tmp_arg) - function(arg)) / dx;
    }

    return result;
}

std::vector<size_t> get_rand_indexes(const size_t length);

std::vector<std::vector<size_t>> get_array_of_rand_indexes(const size_t length, const size_t n);

template <class T>
T factorial(const T& N)
{
    if (N <= 0)
    {
        return 1;
    }
    T result = 1;
    for (T val = 1; val <= N; val++)
    {
        result = result * val;
    }
    return result;
}
}
}
#endif
