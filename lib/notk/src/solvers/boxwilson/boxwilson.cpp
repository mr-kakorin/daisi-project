#include <numeric>

#include <boost/property_tree/ptree.hpp>

#include "../src/tools/tools.h"
#include "boxwilson.h"

#include <notk/notkresult.hpp>

namespace notk
{

class BoxWilsonConfig;

template class BoxWilson<double, double>;
template class BoxWilson<float, float>;

REGISTER_CHILD(BoxWilson<double COMMA double>, IOptimizationStep<double COMMA double>, "Box-Wilson")
REGISTER_CHILD(BoxWilson<float COMMA float>, IOptimizationStep<float COMMA float>, "Box-Wilson")

template <class Targ, class Tfit>
it_res_t<Targ, Tfit> BoxWilson<Targ, Tfit>::do_iteration(const size_t iter_counter,
                                                         const it_res_t<Targ, Tfit>& iter_result,
                                                         const borders_t<Targ>& current_borders)
{
    auto grad = tools::calc_grad(this->m_fitness_function, iter_result.first,
                                 config()->get_gradient_calc_step());

    double h_0 = 1e-5;

    auto apply_grad = [&](const double h) -> std::vector<Targ> {
        std::vector<Targ> cur_x = iter_result.first;
        for (size_t i = 0; i < iter_result.first.size(); i++)
        {
            cur_x[i] = iter_result.first[i] - h * grad[i];
        }
        return cur_x;
    };

    while (this->fitness(apply_grad(h_0)) < iter_result.second)
    {
        h_0 = h_0 * 2.0;
    }

    double h_a = 0;
    double h_b = h_0;

    while (std::abs(h_a - h_b) > h_b * 1e-2)
    {
        if (this->fitness(apply_grad(h_a)) < this->fitness(apply_grad(h_b)))
        {
            h_b = (h_b + h_a) / 2.0;
        }
        else
        {
            h_a = (h_b + h_a) / 2.0;
        }
    }

    return std::make_pair(apply_grad(h_a), this->fitness(apply_grad(h_a)));
}

template <class Targ, class Tfit>
borders_t<Targ> BoxWilson<Targ, Tfit>::squeez_borders(const size_t                iter_counter,
                                                      const it_res_t<Targ, Tfit>& iter_result,
                                                      const borders_t<Targ>& current_borders)
{
    return current_borders;
}

template <class Targ, class Tfit>
bool BoxWilson<Targ, Tfit>::read_config(const boost::property_tree::ptree& config)
{
    return true;
}
template <class Targ, class Tfit>
std::shared_ptr<BoxWilsonConfig> BoxWilson<Targ, Tfit>::config() const
{
    return std::static_pointer_cast<BoxWilsonConfig>(this->m_config);
}
}