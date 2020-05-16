#ifndef NOTK_CONTROLLER_P_H
#define NOTK_CONTROLLER_P_H

#include "states_events.h"

#include <iostream>

#include <boost/msm/back/state_machine.hpp>
#include <boost/msm/front/state_machine_def.hpp>

namespace notk
{

template <class Targ, class Tfit>
class IOptimizationStep;

struct NOTKConfig;

template <class Targ, class Tfit>
class NOTKController_p_
    : public boost::msm::front::state_machine_def<NOTKController_p_<Targ, Tfit>, statemachine::Base>
{

  public:
    typedef statemachine::WaitJSON initial_state;

    NOTKController_p_() noexcept;
    std::shared_ptr<NOTKResults<Targ, Tfit>> result() const noexcept;
    void reset_result() noexcept;

    static void borders_info_msg(const bool correct_borders, const std::string& msg) noexcept;

    bool set_borders(const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
                     const std::vector<Targ>& x_0) noexcept;

    bool set_borders(const std::vector<Targ>& x_left, const std::vector<Targ>& x_right) noexcept;

    bool set_borders(const std::vector<Targ>& x_0, const Targ& delta) noexcept;

    bool set_fitness(const statemachine::set_fitness<Targ, Tfit>& event) noexcept;

    bool set_borders_fitness(const statemachine::set_borders_fitness<Targ, Tfit>& event) noexcept;

    void process(const statemachine::process& event) noexcept;

    typedef NOTKController_p_ p; // makes transition table cleaner

    template <class T>
    void none(const T& val)
    {
    }

    template <typename T1, class Event, typename T2, void (p::*action)(Event const&)>
    using a_row_t = typename boost::msm::front::state_machine_def<
        NOTKController_p_<Targ, Tfit>>::template a_row<T1, Event, T2, action>;

    template <typename T1, class Event, typename T2>
    using _row_t = typename boost::msm::front::state_machine_def<
        NOTKController_p_<Targ, Tfit>>::template _row<T1, Event, T2>;

    template <typename T1, class Event, typename T2, void (p::*action)(Event const&),
              bool (p::*guard)(Event const&)>
    using row_t = typename boost::msm::front::state_machine_def<
        NOTKController_p_<Targ, Tfit>>::template row<T1, Event, T2, action, guard>;

    template <class FSM, class Event>
    void no_transition(Event const& e, FSM& fsm, int state)
    {
        LOG(sev_lvl::error) << "no transition from state: " << fsm.get_state_by_id(state)->name()
                            << " on event: " << e.name();
    }

    struct transition_table
        : boost::mpl::vector<
              //        Start                     Event                   Next                Action
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+
              _row_t<statemachine::WaitJSON, statemachine::json_recieved, statemachine::Configure>,
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+
              _row_t<statemachine::Configure, statemachine::json_non_borders,
                     statemachine::WaitFitnessAndBorders>,
              _row_t<statemachine::Configure, statemachine::json_borders,
                     statemachine::WaitFitness>,
              _row_t<statemachine::Configure, statemachine::invalid_input, statemachine::WaitJSON>,
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+

              _row_t<statemachine::WaitFitness, statemachine::set_fitness<Targ, Tfit>,
                     statemachine::WaitJSON>,
              row_t<statemachine::WaitFitness, statemachine::set_fitness<Targ, Tfit>,
                    statemachine::WaitProcessJSONBorders, &p::none, &p::set_fitness>,
              _row_t<statemachine::WaitFitness, statemachine::json_recieved,
                     statemachine::Configure>,
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+
              a_row_t<statemachine::WaitProcessJSONBorders, statemachine::process,
                      statemachine::WaitFitness, &p::process>,
              _row_t<statemachine::WaitProcessJSONBorders, statemachine::set_fitness<Targ, Tfit>,
                     statemachine::WaitJSON>,
              row_t<statemachine::WaitProcessJSONBorders, statemachine::set_fitness<Targ, Tfit>,
                    statemachine::WaitProcessJSONBorders, &p::none, &p::set_fitness>,
              _row_t<statemachine::WaitProcessJSONBorders, statemachine::json_recieved,
                     statemachine::Configure>,
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+
              _row_t<statemachine::WaitFitnessAndBorders,
                     statemachine::set_borders_fitness<Targ, Tfit>, statemachine::WaitJSON>,
              row_t<statemachine::WaitFitnessAndBorders,
                    statemachine::set_borders_fitness<Targ, Tfit>,
                    statemachine::WaitProcessVectorBorders, &p::none, &p::set_borders_fitness>,
              _row_t<statemachine::WaitFitnessAndBorders, statemachine::set_fitness<Targ, Tfit>,
                     statemachine::WaitJSON>,
              row_t<statemachine::WaitFitnessAndBorders, statemachine::set_fitness<Targ, Tfit>,
                    statemachine::WaitFitnessAndBorders, &p::none, &p::set_fitness>,
              _row_t<statemachine::WaitFitnessAndBorders, statemachine::json_recieved,
                     statemachine::Configure>,
              //  +-------------------------+------------------+------------------------+---------------------+----------------------+
              a_row_t<statemachine::WaitProcessVectorBorders, statemachine::process,
                      statemachine::WaitFitnessAndBorders, &p::process>,
              _row_t<statemachine::WaitProcessVectorBorders, statemachine::set_fitness<Targ, Tfit>,
                     statemachine::WaitJSON>,
              row_t<statemachine::WaitProcessVectorBorders, statemachine::set_fitness<Targ, Tfit>,
                    statemachine::WaitProcessVectorBorders, &p::none, &p::set_fitness>,
              _row_t<statemachine::WaitProcessVectorBorders, statemachine::json_recieved,
                     statemachine::Configure>,
              _row_t<statemachine::WaitProcessVectorBorders,
                     statemachine::set_borders_fitness<Targ, Tfit>, statemachine::WaitJSON>,
              row_t<statemachine::WaitProcessVectorBorders,
                    statemachine::set_borders_fitness<Targ, Tfit>,
                    statemachine::WaitProcessVectorBorders, &p::none, &p::set_borders_fitness>>
    {
    };

    void read_problem_config(const std::shared_ptr<boost::property_tree::ptree>& pt) noexcept;

    commtools::pc_unique_ptr<NOTKConfig>& config();

    const commtools::pc_unique_ptr<NOTKConfig>& config() const;

    void add_observer(const std::shared_ptr<NOTKObserver>& observer);

  protected:
    borders_t<Targ> m_borders;
    std::vector<Targ> m_x0;

    std::shared_ptr<NOTKResults<Targ, Tfit>> m_result;

    commtools::pc_unique_ptr<NOTKConfig> m_config;

    std::list<commtools::pc_unique_ptr<IOptimizationStep<Targ, Tfit>>> m_steps;

    std::list<commtools::pc_shared_ptr<NOTKObserver>> m_observers;

    std::string get_version() noexcept;

    bool read_steps_config(const boost::property_tree::ptree& pt);

    bool set_steps_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function);
};

template <class Targ, class Tfit>
struct NOTKController_p
{
    boost::msm::back::state_machine<NOTKController_p_<Targ, Tfit>> data;
    std::string current_state_name() const
    {
        auto id   = data.current_state()[0];
        auto base = data.get_state_by_id(id);
        return base->name();
    }
};
}

#endif
