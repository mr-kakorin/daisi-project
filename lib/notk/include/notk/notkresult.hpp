#ifndef NOTK_RESULT_H
#define NOTK_RESULT_H

#include <cmath>
#include <cstddef>
#include <vector>

#include <notk/types.hpp>

namespace notk
{
template <class Targ, class Tfit>
class IOptimizationStep;

template <class Targ, class Tfit>
class NOTKController_p_;

template <class Targ, class Tfit>
class NOTKResults
{
    friend class IOptimizationStep<Targ, Tfit>;
    friend class NOTKController_p_<Targ, Tfit>;

  public:
    std::vector<Targ> get_last_argument() const noexcept;

    Tfit get_fitness() const noexcept;

    it_res_t<Targ, Tfit> get_last_it_res() const noexcept;
    it_res_t<Targ, Tfit> get_res(const size_t number) const noexcept;

    template <class Tres>
    std::vector<Tres> get_fit_array() const noexcept;

    std::vector<std::vector<Targ>> get_arg_arr() const noexcept;

    std::vector<std::vector<Targ>> get_arg_non_zero() const noexcept;

    const std::vector<it_res_t<Targ, Tfit>>& result() const noexcept;

  private:
    template <class T1, class T2>
    void add_data(T1&& fitness_argument, T2&& fitnes)
    {
        m_result.emplace_back(std::make_pair(std::forward<std::vector<Targ>>(fitness_argument),
                                             std::forward<Tfit>(fitnes)));
    }

    template <class T>
    void add_data(T&& iter_result)
    {
        m_result.push_back(std::forward<T>(iter_result));
    }

    std::vector<it_res_t<Targ, Tfit>> m_result;
};
}

#endif
