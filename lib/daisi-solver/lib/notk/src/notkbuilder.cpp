#include <notk/notkbuilder.hpp>

namespace notk
{

template <class Targ, class Tfit, typename... Args>
void setter_adapter(std::shared_ptr<NOTKController<Targ, Tfit>>& controller,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                    const Args&... args)
{
    controller->set_borders_fitness(fitness_function, args...);
}

template <class Targ, class Tfit>
void setter_adapter(std::shared_ptr<NOTKController<Targ, Tfit>>& controller,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
{
    controller->set_fitness(fitness_function);
}

template <class Targ, class Tfit, typename... Args>
std::shared_ptr<NOTKController<Targ, Tfit>>
build_notk(const bool is_file, const std::string& config,
           const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
           const Args&... args)
{
    auto NOTKController = std::make_shared<notk::NOTKController<Targ, Tfit>>();

    if (is_file)
    {
        NOTKController->set_problem_config(config);
    }
    else
    {
        NOTKController->set_problem_config_str(config);
    }

    setter_adapter(NOTKController, fitness_function, args...);

    return NOTKController->ready_to_process() ? NOTKController : nullptr;
}

template struct NOTKBuilder<double, double>;
template struct NOTKBuilder<float, float>;

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_file_config(
    const std::string&                                   config_file_path,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right, const std::vector<Targ>& x_0)
{

    return build_notk(true, config_file_path, fitness_function, x_left, x_right, x_0);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_file_config(
    const std::string&                                   config_file_path,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right)
{
    return build_notk(true, config_file_path, fitness_function, x_left, x_right);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_file_config(
    const std::string&                                   config_file_path,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_0, const Targ& delta)
{
    return build_notk(true, config_file_path, fitness_function, x_0, delta);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_file_config(
    const std::string&                                   config_file_path,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
{
    return build_notk(true, config_file_path, fitness_function);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_string_config(
    const std::string&                                   config,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right, const std::vector<Targ>& x_0)
{
    return build_notk(false, config, fitness_function, x_left, x_right, x_0);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_string_config(
    const std::string&                                   config,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right)
{
    return build_notk(false, config, fitness_function, x_left, x_right);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_string_config(
    const std::string&                                   config,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
    const std::vector<Targ>& x_0, const Targ& delta)
{
    return build_notk(false, config, fitness_function, x_0, delta);
}

template <class Targ, class Tfit>
std::shared_ptr<NOTKController<Targ, Tfit>> NOTKBuilder<Targ, Tfit>::build_notk_string_config(
    const std::string&                                   config,
    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function)
{
    return build_notk(false, config, fitness_function);
}
}
