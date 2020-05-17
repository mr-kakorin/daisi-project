#include <numeric>

#include <serreflection/read_json.hpp>

#include "../src/tools/tools.h"
#include "heavyball.h"
#include "heavyballconfig.h"

#include <notk/notkresult.hpp>

namespace notk
{

template class HeavyBall<double, double>;
template class HeavyBall<float, float>;

REGISTER_CHILD(HeavyBall<double COMMA double>, IOptimizationStep<double COMMA double>, "Heavy ball")
REGISTER_CHILD(HeavyBall<float COMMA float>, IOptimizationStep<float COMMA float>, "Heavy ball")

template <class Targ, class Tfit>
it_res_t<Targ, Tfit> HeavyBall<Targ, Tfit>::do_iteration(const size_t                iter_counter,
                                                         const it_res_t<Targ, Tfit>& iter_result,
                                                         const borders_t<Targ>& current_borders)
{
    auto grad = tools::calc_grad(this->m_fitness_function, iter_result.first,
                                 config()->get_gradient_calc_step());

    double h_0 = config()->get_gradient_calc_step();

    auto start_weight = config()->get_weight();

    auto apply_grad = [&](const double h, bool& border) -> std::vector<Targ> {
        std::vector<Targ> cur_x = iter_result.first;
        for (size_t i = 0; i < iter_result.first.size(); i++)
        {
            Targ ball_coef = 0;
            if (last_x.size())
            {
                ball_coef = start_weight * (iter_result.first[i] - last_x[i]);
            }
            cur_x[i] = iter_result.first[i] - h * grad[i] + ball_coef;
            if (cur_x[i] > current_borders.second[i])
            {
                border   = false;
                cur_x[i] = current_borders.second[i];
            }
            if (cur_x[i] < current_borders.first[i])
            {
                border   = false;
                cur_x[i] = current_borders.first[i];
            }
        }
        return cur_x;
    };

    bool border = true;

    bool   flag = true;
    double h_a;
    double h_b;

    while (flag)
    {
        auto x2 = this->fitness(apply_grad(2 * h_0, border));

        if (x2 < iter_result.second)
        {
            while (this->fitness(apply_grad(h_0, border)) < iter_result.second && border)
            {
                h_0 = h_0 * 2.0;
            }
            h_a  = h_0 / 2.0;
            h_b  = h_0;
            flag = false;
        }
        else
        {
            bool ok = true;
            while (this->fitness(apply_grad(h_0, border)) > iter_result.second && border)
            {
                if (h_0 > -1e-200 || start_weight == 0)
                {
                    h_0 = h_0 / 2.0;
                }
                else
                {
                    start_weight = start_weight / 2.0;

                    h_0 = config()->get_gradient_calc_step();
                    ok  = false;
                    break;
                }
            }
            if (ok)
            {
                flag = false;
                h_a  = 0;
                h_b  = 2 * h_0;
            }
        }
    }

    auto local_fitness = [&](const double h) -> double {
        bool border = true;
        return this->fitness(apply_grad(h, border));
    };

    m_searcher1d->set_fitness(local_fitness);

    borders_t<Targ>   borders_1d    = std::make_pair(std::vector<Targ>{static_cast<Targ>(h_a)},
                                                std::vector<Targ>{static_cast<Targ>(h_b)});
    std::vector<Targ> x_current_tmp = {static_cast<Targ>((h_a + h_b) / 2.0)};

    bool flag_abort = true;

    NOTKResults<Targ, Tfit> result_tmp;

    m_searcher1d->opt_routine(result_tmp, x_current_tmp, borders_1d, 1, true, flag_abort, false);

    auto x_best = apply_grad(result_tmp.get_last_argument()[0], border);

    last_x = iter_result.first;

    return std::make_pair(x_best, this->fitness(x_best));
}

template <class Targ, class Tfit>
borders_t<Targ> HeavyBall<Targ, Tfit>::squeez_borders(const size_t                iter_counter,
                                                      const it_res_t<Targ, Tfit>& iter_result,
                                                      const borders_t<Targ>&      current_borders)
{
    return current_borders;
}
template <class Targ, class Tfit>
bool HeavyBall<Targ, Tfit>::read_config(const boost::property_tree::ptree& config)
{
    this->m_config = srfl::read_json<HeavyBallConfig>(config);

    try
    {
        auto pt_searcher1d = config.get_child("searcher1d");

        m_searcher1d = IOptimizationStep<Targ, Tfit>::create_step(pt_searcher1d);

        m_searcher1d->read_config(pt_searcher1d);

        return this->m_config && m_searcher1d->config() ? true : false;
    }
    catch (const std::exception& ex)
    {
        BL_ERROR() << "Error read searcher 1d config: " << ex.what();
        return false;
    }
}

template <class Targ, class Tfit>
std::shared_ptr<HeavyBallConfig> HeavyBall<Targ, Tfit>::config() const
{
    return std::static_pointer_cast<HeavyBallConfig>(this->m_config);
}
} // namespace notk