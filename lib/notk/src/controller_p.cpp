#include <utility>

#include <common_tools/boostlog.hpp>
#include <common_tools/json_helper.h>
#include <serreflection/read_json.hpp>

#include <notk/controller.hpp>

#include "base/ioptistep.h"
#include "config.h"
#include "controller_p.h"
#include "tools/config_parser_helper.h"

namespace notk
{

const std::string error_msg_controller_prefix = "Read config of optimization::error: ";

template class NOTKController_p_<double, double>;
template class NOTKController_p_<float, float>;

template <class Targ, class Tfit>
NOTKController_p_<Targ, Tfit>::NOTKController_p_() noexcept : m_result(nullptr), m_config(nullptr)
{
}

template <class Targ, class Tfit>
std::string NOTKController_p_<Targ, Tfit>::get_version() noexcept
{
    return "0.3";
}

template <class Targ, class Tfit>
commtools::pc_unique_ptr<NOTKConfig>& NOTKController_p_<Targ, Tfit>::config()
{
    return m_config;
}

template <class Targ, class Tfit>
const commtools::pc_unique_ptr<NOTKConfig>& NOTKController_p_<Targ, Tfit>::config() const
{
    return m_config;
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKResults<Targ, Tfit>> NOTKController_p_<Targ, Tfit>::result() const
    noexcept
{
    return m_result;
}

template <class Targ, class Tfit>
void NOTKController_p_<Targ, Tfit>::reset_result() noexcept
{
    m_result = nullptr;
}

template <class Targ, class Tfit>
void NOTKController_p_<Targ, Tfit>::process(const statemachine::process& event) noexcept
{
    reset_result();
    try
    {
        borders_t<Targ> borders(m_borders);
        std::vector<Targ> x0(m_x0);

        BL_INFO() << "start calculations";

        m_result = std::make_shared<NOTKResults<Targ, Tfit>>();

        for (auto& step : m_steps)
        {
            step->opt_routine(*m_result, x0, borders, m_config->log_each_iteration,
                              m_config->allow_maximization, event.m_flag_abort, true);
        }
        if (m_result->m_result.empty())
        {
            BL_ERROR() << "no solution found";
            reset_result();
            return;
        }
        BL_INFO() << "calculations is done, best fitness = " << m_result->get_fitness();
    }
    catch (const std::exception& ex)
    {
        BL_ERROR() << "calculation error: " << ex.what();
        reset_result();
    }
    catch (...)
    {
        BL_ERROR() << "unknown calculation error";
        reset_result();
    }
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::read_steps_config(const boost::property_tree::ptree& pt)
{

    auto error_msg_prefix = "Error read solver config";
    return commtools::read_property_tree(
        pt, error_msg_prefix, [&](const boost::property_tree::ptree& pt) -> bool {
            m_steps.clear();

            bool result = true;

            size_t n_steps = 0;
            for (auto& item : pt.get_child("methods"))
            {
                boost::property_tree::ptree pt = item.second;

                auto step = IOptimizationStep<Targ, Tfit>::create_step(pt);

                if (step == nullptr)
                {
                    result = false;
                }
                else
                {
                    if (!step->read_config(pt))
                    {
                        result = false;
                    }
                    m_steps.push_back(std::move(step));
                    n_steps++;
                }
            }
            if (result == true && n_steps == 0)
            {
                BL_ERROR() << error_msg_prefix << "no steps";
                result = false;
            }
            if (result == true && n_steps != 0)
            {
                BL_INFO() << "configuration for " << n_steps << " steps success";
            }
            return result;
        });
}

template <class Targ, class Tfit>
void NOTKController_p_<Targ, Tfit>::read_problem_config(
    const std::shared_ptr<boost::property_tree::ptree>& pt) noexcept
{
    BL_INFO() << "set NOTK configuration";

    const std::string err_msg = "set NOTK configuration fail";

    m_config = nullptr;

    if (!pt)
    {
        BL_ERROR() << err_msg;
        return;
    }

    m_config = srfl::read_json<NOTKConfig>(*pt);

    if (!m_config)
    {
        BL_ERROR() << err_msg;
    }
    else if (notk::InitialBordersType::MANUAL == m_config->borders_type)
    {
        try
        {
            auto pt_borders = pt->get_child("borders");

            auto borders = srfl::read_json<NOTKBorders>(pt_borders);

            if (!borders)
            {
                throw std::runtime_error("incorrect borders config");
            }

            auto resize_borders = [&](std::vector<double>& border) -> std::vector<Targ> {

                std::vector<Targ> result(borders->problem_dimension);

                int cp_len = std::min(borders->problem_dimension, static_cast<int>(border.size()));

                std::copy(border.begin(), border.begin() + cp_len, result.begin());

                if (int(border.size()) < borders->problem_dimension)
                {
                    auto fil_vall = static_cast<Targ>(border.back());
                    std::fill(result.begin() + border.size(), result.end(), fil_vall);
                    BL_INFO() << "Atomatical border filling, value = " << fil_vall;
                }
                return result;
            };

            std::vector<Targ> x_left  = resize_borders(borders->x_left);
            std::vector<Targ> x_right = resize_borders(borders->x_right);
            std::vector<Targ> x_0     = resize_borders(borders->x_0);

            if (!set_borders(x_left, x_right, x_0))
            {
                m_config = nullptr;
            }
        }
        catch (const std::exception& ex)
        {
            BL_ERROR() << "Error read optimization borders " << ex.what();
            m_config = nullptr;
        }
    }
    else if (notk::InitialBordersType::VECTOR == m_config->borders_type)
    {
        BL_INFO() << "wait set optimization borders from vectors";
    }
    else
    {
        BL_ERROR() << "unexpected borders type";
        m_config = nullptr;
    }

    if (m_config)
    {
        if (!read_steps_config(*pt))
        {
            m_config = nullptr;
        }
    }
}

template <class Targ, class Tfit>
void NOTKController_p_<Targ, Tfit>::borders_info_msg(const bool correct_borders,
                                                     const std::string& msg) noexcept
{
    if (correct_borders)
    {
        BL_INFO() << "set optimization borders success";
    }
    else
    {
        BL_INFO() << "set optimization borders fail: " << msg;
    }
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_steps_fitness(
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
{
    bool result = true;
    for (auto& step : m_steps)
    {
        if (!step->set_fitness(fitness_function))
        {
            result = false;
        }
    }
    return result;
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_borders(const std::vector<Targ>& x_left,
                                                const std::vector<Targ>& x_right,
                                                const std::vector<Targ>& x_0) noexcept
{
    BL_INFO() << "set optimization borders from x_left, x_right and x_0";

    bool correct_borders = true;

    std::string msg = "";

    CHECKVALUE(x_left.size() != x_0.size(), correct_borders, msg);
    CHECKVALUE(x_right.size() != x_0.size(), correct_borders, msg);

    if (!correct_borders)
    {
        borders_info_msg(correct_borders, msg);

        return correct_borders;
    }

    for (size_t i = 0; i < x_0.size(); i++)
    {
        CHECKVALUE(x_left[i] > x_0[i] || x_0[i] > x_right[i], correct_borders, msg);
    }
    borders_info_msg(correct_borders, msg);

    if (!correct_borders)
    {
        return correct_borders;
    }

    m_borders = std::make_pair(x_left, x_right);
    m_x0      = x_0;

    return correct_borders;
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_borders(const std::vector<Targ>& x_left,
                                                const std::vector<Targ>& x_right) noexcept
{
    bool correct_borders = true;

    std::string msg = "";

    CHECKVALUE(x_left.size() != x_right.size(), correct_borders, msg);

    if (!correct_borders)
    {
        borders_info_msg(correct_borders, msg);
        return correct_borders;
    }

    for (size_t i = 0; i < x_left.size(); i++)
    {
        CHECKVALUE(x_left[i] >= x_right[i], correct_borders, msg);
    }

    borders_info_msg(correct_borders, msg);

    if (!correct_borders)
    {
        return correct_borders;
    }

    auto problem_dim = x_left.size();
    m_borders        = std::make_pair(x_left, x_right);

    BL_WARNING() << "set_borders::force set x_0 as average";
    m_x0.resize(problem_dim);
    for (size_t i = 0; i < problem_dim; i++)
    {
        m_x0[i] = (x_left[i] + x_right[i]) / 2.0;
    }
    return correct_borders;
}
template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_borders(const std::vector<Targ>& x_0,
                                                const Targ& delta) noexcept
{
    BL_INFO() << "set optimization borders from x0 and delta";

    bool correct_borders = true;

    std::string msg = "";

    CHECKVALUE(delta <= 0, correct_borders, msg);

    auto problem_dim = static_cast<int>(x_0.size());

    borders_info_msg(correct_borders, msg);

    if (!correct_borders)
    {
        return correct_borders;
    }

    m_borders.first.resize(problem_dim);
    m_borders.second.resize(problem_dim);
    for (int i = 0; i < problem_dim; i++)
    {
        m_borders.first[i]  = x_0[i] - delta;
        m_borders.second[i] = x_0[i] + delta;
    }
    return correct_borders;
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_borders_fitness(
    const statemachine::set_borders_fitness<Targ, Tfit>& event) noexcept
{
    reset_result();

    if (!set_fitness(statemachine::set_fitness<Targ, Tfit>(event.m_fitness_function)))
    {
        return false;
    }
    switch (event.m_type)
    {
    case statemachine::set_borders_fitness<Targ, Tfit>::InputBordersType::XL_XR:
        return set_borders(event.m_x_left, event.m_x_right);
    case statemachine::set_borders_fitness<Targ, Tfit>::InputBordersType::XL_XR_X0:
        return set_borders(event.m_x_left, event.m_x_right, event.m_x_0);
    case statemachine::set_borders_fitness<Targ, Tfit>::InputBordersType::X0_DELTA:
        return set_borders(event.m_x_0, event.m_delta);
    default:
        BL_ERROR() << "unexpected borders type";
        return false;
    }
}

template <class Targ, class Tfit>
bool NOTKController_p_<Targ, Tfit>::set_fitness(
    const statemachine::set_fitness<Targ, Tfit>& event) noexcept
{
    reset_result();

    bool correct_fitness = true;

    BL_INFO() << "set optimization fitness";

    if (!set_steps_fitness(event.m_fitness_function))
    {
        correct_fitness = false;
    }
    if (correct_fitness)
    {
        BL_INFO() << "set optimization fitness success";
    }
    else
    {
        BL_INFO() << "set optimization fitness fail";
    }
    return correct_fitness;
}

template <class Targ, class Tfit>
void NOTKController_p_<Targ, Tfit>::add_observer(const std::shared_ptr<NOTKObserver>& observer)
{
    m_observers.push_back(observer);
}

}
