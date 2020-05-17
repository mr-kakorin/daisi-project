#include <common_tools/boostlog.hpp>
#include <common_tools/json_helper.h>

#include <notk/notkresult.hpp>
#include <notk/types.hpp>

#include "iconfig.h"
#include "ioptistep.h"

namespace notk
{
template class IOptimizationStep<double, double>;
template class IOptimizationStep<float, float>;

template <class Targ, class Tfit>
IOptimizationStep<Targ, Tfit>::~IOptimizationStep() = default;

template <class Targ, class Tfit>
IOptimizationStep<Targ, Tfit>::IOptimizationStep() : m_config(nullptr)
{
}

template <class Targ, class Tfit>
std::unique_ptr<IOptimizationStep<Targ, Tfit>>
IOptimizationStep<Targ, Tfit>::create_step(const boost::property_tree::ptree& pt)
{
    std::string type = pt.get<std::string>("type", "");

    auto result = make_unique_child(type);

    if (!result)
    {
        BL_ERROR() << "unexpected solver type: " << type;
    }

    return result;
}

template <class Targ, class Tfit>
std::shared_ptr<BaseOptConfig> IOptimizationStep<Targ, Tfit>::config() const
{
    return m_config;
}

template <class Targ, class Tfit>
bool IOptimizationStep<Targ, Tfit>::set_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
{
    isMultidim         = true;
    m_fitness_function = fitness_function;
    return true;
}

template <class Targ, class Tfit>
bool IOptimizationStep<Targ, Tfit>::set_fitness(
    const std::function<Tfit(const Targ&)>& fitness_function)
{
    m_fitness_function1d = fitness_function;
    isMultidim           = false;
    return true;
}

template <class Targ, class Tfit>
void IOptimizationStep<Targ, Tfit>::opt_routine(
    NOTKResults<Targ, Tfit>& result, std::vector<Targ>& x_best, borders_t<Targ>& current_borders,
    const unsigned log_each_iteration, const bool allow_maximization, bool& flag_abort, bool do_log,
    const std::list<commtools::pc_shared_ptr<NOTKObserver>>& observers)
{
    check_data();

    Targ fit_best = fitness(x_best);

    size_t iter_counter = 0;

    if (do_log)
    {
        BL_TRACE_CH("NOTK_RES") << fit_best;
        BL_TRACE() << "Iteration: " << iter_counter << ", fitness: " << fit_best;
    }

    result.add_data(x_best, fit_best);

    result.add_data(do_preprocess(result.get_last_it_res(), current_borders));

    if (do_log)
    {
        BL_TRACE_CH("NOTK_RES") << fit_best;
        BL_TRACE() << "Preprocess, fitness: " << result.get_fitness();
    }

    bool accuracy_achieved =
        m_config->check_accuracy(result.get_res(0), result.get_res(1), false, true);

    while (((!accuracy_achieved) && iter_counter < m_config->get_maxNumIterations()) && flag_abort)
    {
        iter_counter++;

        it_res_t<Targ, Tfit> iter_result =
            do_iteration(iter_counter, result.get_last_it_res(), current_borders);

        if (!allow_maximization)
        {
            if (iter_result.second > result.get_fitness())
            {
                iter_result = result.get_last_it_res();
            }
        }
        bool log_local = iter_counter % log_each_iteration == 0;

        accuracy_achieved = m_config->check_accuracy(result.get_last_it_res(), iter_result,
                                                     do_log && log_local, false);

        result.add_data(iter_result);

        current_borders = squeez_borders(iter_counter, result.get_last_it_res(), current_borders);

        x_best = result.get_last_argument();

        if (do_log)
        {
            BL_TRACE_CH("NOTK_RES") << iter_result.second;

            std::stringstream log_msg;

            if (!observers.empty() || log_local)
            {
                log_msg << "Iteration: " << iter_counter << ", fitness: " << result.get_fitness();
            }

            for (const auto& obs : observers)
            {
                obs->handle_event(log_msg.str(),
                                  iter_counter / double(m_config->get_maxNumIterations()));
            }

            if (log_local)
            {
                BL_TRACE() << log_msg.str();
            }
        }
    }
}

template <class Targ, class Tfit>
Tfit IOptimizationStep<Targ, Tfit>::fitness(const std::vector<Targ>& arg) const
{
    if (0 == arg.size())
    {
        throw std::runtime_error("IOptimizationStep::zero size of input argument");
    }
    return isMultidim ? m_fitness_function(arg) : m_fitness_function1d(arg[0]);
}

template <class Targ, class Tfit>
void IOptimizationStep<Targ, Tfit>::check_data() const
{
}

template <class Targ, class Tfit>
it_res_t<Targ, Tfit>
IOptimizationStep<Targ, Tfit>::do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                             const borders_t<Targ>& current_borders)
{
    return iter_result;
}
}
