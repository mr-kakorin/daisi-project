#ifndef NOTK_STATES_EVENTS_H
#define NOTK_STATES_EVENTS_H

#include <boost/msm/back/state_machine.hpp>
#include <boost/msm/front/state_machine_def.hpp>
#include <boost/property_tree/ptree.hpp>

#include <common_tools/boostlog.hpp>

#include "config.h"
#include "enums.h"

namespace notk
{
namespace statemachine
{

namespace statenames
{
const static std::string WaitJSON                 = "WaitJSON";
const static std::string Configure                = "Configure";
const static std::string WaitFitnessAndBorders    = "WaitFitnessAndBorders";
const static std::string WaitFitness              = "WaitFitness";
const static std::string WaitProcessJSONBorders   = "WaitProcessJSONBorders";
const static std::string WaitProcessVectorBorders = "WaitProcessVectorBorders";
}

// events
struct BaseEvent
{
    virtual std::string name() const = 0;
};

struct json_recieved : BaseEvent
{
    json_recieved(const std::shared_ptr<boost::property_tree::ptree>& json) : m_json(json)
    {
    }
    std::shared_ptr<boost::property_tree::ptree> m_json;

    std::string name() const
    {
        return "json_recieved";
    }
};

template <class Tfit, class Targ>
struct set_fitness : BaseEvent
{
    set_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
        : m_fitness_function(fitness_function)
    {
    }
    const std::function<Tfit(const std::vector<Targ>&)>& m_fitness_function;

    std::string name() const
    {
        return "set_fitness";
    }
};

template <class Tfit, class Targ>
struct set_borders_fitness : BaseEvent
{
    std::string name() const
    {
        return "set_borders_fitness";
    }

    enum class InputBordersType
    {
        XL_XR,
        XL_XR_X0,
        X0_DELTA
    };

    set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                        const std::vector<Targ>& x_left, const std::vector<Targ>& x_right)
        : m_fitness_function(fitness_function), m_x_left(x_left), m_x_right(x_right), m_x_0(x_left),
          m_type(InputBordersType::XL_XR)
    {
    }

    set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                        const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
                        const std::vector<Targ>& x_0)
        : m_fitness_function(fitness_function), m_x_left(x_left), m_x_right(x_right), m_x_0(x_0),
          m_type(InputBordersType::XL_XR_X0)
    {
    }

    set_borders_fitness(const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                        const std::vector<Targ>& x_0, const double delta)
        : m_fitness_function(fitness_function), m_x_left(x_0), m_x_right(x_0), m_x_0(x_0),
          m_delta(delta), m_type(InputBordersType::X0_DELTA)
    {
    }

    const std::function<Tfit(const std::vector<Targ>&)>& m_fitness_function;

    const std::vector<Targ>& m_x_left;
    const std::vector<Targ>& m_x_right;
    const std::vector<Targ>& m_x_0;

    double m_delta;

    InputBordersType m_type;
};
struct json_non_borders : BaseEvent
{
    std::string name() const
    {
        return "json_non_borders";
    }
};
struct json_borders : BaseEvent
{
    std::string name() const
    {
        return "json_borders";
    }
};
struct invalid_input : BaseEvent
{
    std::string name() const
    {
        return "invalid_input";
    }
};
struct process : BaseEvent
{
    process(bool& flag_abort) : m_flag_abort(flag_abort)
    {
    }
    bool& m_flag_abort;

    std::string name() const
    {
        return "process";
    }
};
// states
struct Base
{
    virtual std::string name() const
    {
        return "Base";
    }
};

struct WaitJSON : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        BL_INFO() << "enter state: " << name();
    }
    virtual std::string name() const override final
    {
        return statenames::WaitJSON;
    }
};

struct Configure : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        fsm.reset_result();

        BL_INFO() << "enter state: " << name();

        fsm.read_problem_config(ev.m_json);

        if (!fsm.config().get())
        {
            fsm.process_event(invalid_input());
        }
        else if (fsm.config()->borders_type == notk::InitialBordersType::MANUAL)
        {
            fsm.process_event(json_borders());
        }
        else if (fsm.config()->borders_type == notk::InitialBordersType::VECTOR)
        {
            fsm.process_event(json_non_borders());
        }
        else
        {
            throw std::runtime_error("Configure state: Unexpected borders type");
        }
    }
    virtual std::string name() const override final
    {
        return statenames::Configure;
    }
};

struct WaitFitnessAndBorders : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        BL_INFO() << "enter state: " << name();
    }

    virtual std::string name() const override final
    {
        return statenames::WaitFitnessAndBorders;
    }
};

// The list of FSM states
struct WaitFitness : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        BL_INFO() << "enter state: " << name();
    }
    virtual std::string name() const override final
    {
        return statenames::WaitFitness;
    }
};

struct WaitProcessJSONBorders : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        BL_INFO() << "enter state: " << name();
    }
    virtual std::string name() const override final
    {
        return statenames::WaitProcessJSONBorders;
    }
};

struct WaitProcessVectorBorders : public boost::msm::front::state<Base>
{
    template <class Event, class FSM>
    void on_entry(const Event& ev, FSM& fsm)
    {
        BL_INFO() << "enter state: " << name();
    }
    virtual std::string name() const override final
    {
        return statenames::WaitProcessVectorBorders;
    }
};
}
}

#endif
