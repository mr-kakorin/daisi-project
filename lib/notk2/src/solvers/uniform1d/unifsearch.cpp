#include <serreflection/read_json.hpp>

#include "unifsearch.h"

#include "unifsearchconfig.h"

namespace notk
{
template class UniformSearch1d<double, double>;
template class UniformSearch1d<float, float>;

REGISTER_CHILD(UniformSearch1d<double COMMA double>, IOptimizationStep<double COMMA double>)
REGISTER_CHILD(UniformSearch1d<float COMMA float>, IOptimizationStep<float COMMA float>)

template <class Targ, class Tfit>
it_res_t<Targ, Tfit>
UniformSearch1d<Targ, Tfit>::do_iteration(const size_t                 iter_counter,
                                          const it_res_t<Targ, Tfit>&  iter_result,
                                          const borders_t<Targ>& current_borders)
{
    size_t n_divisions_iter;

    Targ dx = calc_dx(iter_counter, current_borders, n_divisions_iter);

    Tfit fit_best = std::numeric_limits<Tfit>::infinity();

    Targ x_best;

    for (size_t i = 0; i < n_divisions_iter; i++)
    {
        Targ x   = current_borders.first[0] + i * dx;
        Tfit fit = this->fitness(std::vector<Targ>{x});

        if (fit < fit_best)
        {
            fit_best = fit;
            x_best   = x;
        }
    }

    it_res_t<Targ, Tfit> result = std::make_pair(std::vector<Targ>{x_best}, fit_best);
    return result;
}

template <class Targ, class Tfit>
borders_t<Targ>
UniformSearch1d<Targ, Tfit>::squeez_borders(const size_t                 iter_counter,
                                            const it_res_t<Targ, Tfit>&  iter_result,
                                            const borders_t<Targ>& current_borders)
{
    size_t n_divisions_iter;
    Targ   dx = calc_dx(iter_counter, current_borders, n_divisions_iter);
    Targ   x0 = iter_result.first[0];
    borders_t<Targ> result =
        std::make_pair(std::vector<Targ>{std::max(x0 - dx, current_borders.first[0])},
                       std::vector<Targ>{std::min(x0 + dx, current_borders.second[0])});

    return result;
}

template <class Targ, class Tfit>
Targ UniformSearch1d<Targ, Tfit>::calc_dx(const size_t                 iter_counter,
                                          const borders_t<Targ>& current_borders,
                                          size_t& n_divisions_iter) const
{
    n_divisions_iter = config()->get_n_divisions_first();

    if (iter_counter > 1)
    {
        n_divisions_iter = config()->get_n_divisions();
    }

    Targ dx = (current_borders.second[0] - current_borders.first[0]) /
              (static_cast<Targ>(n_divisions_iter) - 1.0);

    return dx;
}
template <class Targ, class Tfit>
bool UniformSearch1d<Targ, Tfit>::read_config(const boost::property_tree::ptree& config)
{
    this->m_config = srfl::read_json<UniformSearch1dConfig>(config);

    return this->m_config ? true : false;
}
template <class Targ, class Tfit>
std::shared_ptr<UniformSearch1dConfig> UniformSearch1d<Targ, Tfit>::config() const
{
    return std::static_pointer_cast<UniformSearch1dConfig>(this->m_config);
}
}