#include <boostlog0/boostlog.h>
#include <common_tools/json_helper.h>

#include <notk/controller.hpp>

#include "base/ioptistep.h"
#include "controller_p.h"
#include "states_events.h"
#include "tools/config_parser_helper.h"

namespace notk
{
bool interpret_transition_result(const std::string& name)
{
    return name != statemachine::statenames::WaitJSON;
}

bool interpret_transition_result(boost::msm::back::HandledEnum transition_result)
{
    return transition_result == boost::msm::back::HandledEnum::HANDLED_TRUE;
}

template <class Tfsm, class Teven>
bool process_and_check_event(Tfsm& fsm, const Teven& event)
{
    if (!interpret_transition_result(fsm->data.process_event(event)))
    {
        return false;
    }
    return interpret_transition_result(fsm->current_state_name());
}

template class NOTKController<double, double>;
template class NOTKController<float, float>;

template <class Targ, class Tfit>
NOTKController<Targ, Tfit>::~NOTKController() = default;

template <class Targ, class Tfit>
NOTKController<Targ, Tfit>::NOTKController() noexcept
    : m_ipml(std::make_unique<NOTKController_p<Targ, Tfit>>())
{
    m_ipml->data.start();
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKResults<Targ, Tfit>>
NOTKController<Targ, Tfit>::process(bool& flag_abort) noexcept
{
    if (!interpret_transition_result(m_ipml->data.process_event(statemachine::process(flag_abort))))
    {
        return nullptr;
    }

    return m_ipml->data.result();
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_problem_config(
    const std::shared_ptr<boost::property_tree::ptree>& json) noexcept
{
    return process_and_check_event(m_ipml, statemachine::json_recieved(json));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_problem_config(const std::string& path) noexcept
{
    return set_problem_config(commtools::readJSONFile(path));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_problem_config_str(const std::string& json_string) noexcept
{
    return set_problem_config(commtools::readJSONString(json_string));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::ready_to_process() const noexcept
{
    return m_ipml->current_state_name() == statemachine::statenames::WaitProcessJSONBorders ||
           m_ipml->current_state_name() == statemachine::statenames::WaitProcessVectorBorders;
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_borders_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
    const std::vector<Targ>& x_0) noexcept
{
    return process_and_check_event(m_ipml, statemachine::set_borders_fitness<Targ, Tfit>(
                                               fitness_function, x_left, x_right, x_0));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_borders_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right) noexcept
{
    return process_and_check_event(
        m_ipml, statemachine::set_borders_fitness<Targ, Tfit>(fitness_function, x_left, x_right));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_borders_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_0, const Targ& delta) noexcept
{
    return process_and_check_event(
        m_ipml, statemachine::set_borders_fitness<Targ, Tfit>(fitness_function, x_0, delta));
}

template <class Targ, class Tfit>
bool NOTKController<Targ, Tfit>::set_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function) noexcept
{
    return process_and_check_event(m_ipml, statemachine::set_fitness<Targ, Tfit>(fitness_function));
}

template <class Targ, class Tfit>
void NOTKController<Targ, Tfit>::add_observer(const std::shared_ptr<NOTKObserver>& observer)
{
    m_ipml->data.add_observer(observer);
}

}
