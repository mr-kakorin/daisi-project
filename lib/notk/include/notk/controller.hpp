#ifndef NOTK_CONTROLLER_H
#define NOTK_CONTROLLER_H

#include <boost/thread/future.hpp>

#include <functional>
#include <memory>

#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>
#include <boost/thread/future.hpp>

#include <common_tools/propagate_const.hpp>

#include <notk/notkresult.hpp>
#include <notk/types.hpp>
#include <notk/observer.hpp>

namespace notk
{

template <class Targ, class Tfit>
class NOTKController_p;

template <class Targ, class Tfit>

class NOTKController
{

  public:
    virtual ~NOTKController();
    NOTKController() noexcept;

    std::shared_ptr<NOTKResults<Targ, Tfit>> process(bool& flag_abort) noexcept;

    bool set_problem_config(const std::string& path) noexcept;
    bool set_problem_config_str(const std::string& json_string) noexcept;
    bool set_problem_config(const std::shared_ptr<boost::property_tree::ptree>& json) noexcept;

    bool ready_to_process() const noexcept;

    bool set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                             const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
                             const std::vector<Targ>& x_0) noexcept;

    bool set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                             const std::vector<Targ>&                             x_left,
                             const std::vector<Targ>&                             x_right) noexcept;

    bool set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                             const std::vector<Targ>& x_0, const Targ& delta) noexcept;

    bool
    set_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function) noexcept;

    void add_observer(const std::shared_ptr<NOTKObserver>& observer);

  private:
    commtools::pc_unique_ptr<NOTKController_p<Targ, Tfit>> m_ipml;
};
}

#endif
