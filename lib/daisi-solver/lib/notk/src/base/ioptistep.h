#ifndef NOTK_I_OPTIMIZATIONSTEP_H
#define NOTK_I_OPTIMIZATIONSTEP_H

#include <functional>
#include <memory>
#include <vector>
#include <list>

#include <boost/property_tree/ptree.hpp>

#include <common_tools/child_factory.hpp>
#include <common_tools/propagate_const.hpp>

#include <notk/observer.hpp>
#include <notk/types.hpp>

namespace notk
{
template <class Targ, class Tfit>
class NOTKResults;
class BaseOptConfig;

template <class Targ, class Tfit>
class IOptimizationStep
{
  public:
    virtual ~IOptimizationStep();

    BUILD_CHILD_FACTORY(std::string, IOptimizationStep<Targ COMMA Tfit>)

    IOptimizationStep();

    virtual bool read_config(const boost::property_tree::ptree& config) = 0;

    bool set_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function);

    bool set_fitness(const std::function<Tfit(const Targ&)>& fitness_function);

    void opt_routine(NOTKResults<Targ, Tfit>& result, std::vector<Targ>& x_best,
                     borders_t<Targ>& current_borders, const unsigned log_each_iteration,
                     const bool allow_maximization, bool& flag_abort, bool do_log = true,
                     const std::list<commtools::pc_shared_ptr<NOTKObserver>>& observers = {});

    std::shared_ptr<BaseOptConfig> config() const;

    static std::unique_ptr<IOptimizationStep<Targ, Tfit>>
    create_step(const boost::property_tree::ptree& pt);

  protected:
    Tfit fitness(const std::vector<Targ>& arg) const;

    virtual void check_data() const;

    virtual it_res_t<Targ, Tfit> do_iteration(const size_t                iter_counter,
                                              const it_res_t<Targ, Tfit>& iter_result,
                                              const borders_t<Targ>& current_borders) = 0;

    virtual it_res_t<Targ, Tfit> do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                               const borders_t<Targ>& current_borders);

    virtual borders_t<Targ> squeez_borders(const size_t iter_counter,
                                           const it_res_t<Targ, Tfit>& iter_result,
                                           const borders_t<Targ>& current_borders) = 0;

    std::shared_ptr<BaseOptConfig> m_config;

    bool isMultidim;

    borders_t<Targ> m_borders;

    std::function<Tfit(const std::vector<Targ>&)> m_fitness_function;

    std::function<Tfit(const Targ&)> m_fitness_function1d;
};
}

#endif
