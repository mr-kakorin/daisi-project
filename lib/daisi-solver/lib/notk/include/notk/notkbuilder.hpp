#ifndef NOTK_BUILDER_HPP
#define NOTK_BUILDER_HPP

#include <notk/controller.hpp>

namespace notk
{
template <class Targ, class Tfit>
struct NOTKBuilder
{
    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_file_config(const std::string&                                   config_file_path,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
                    const std::vector<Targ>& x_0);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_file_config(const std::string&                                   config_file_path,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                    const std::vector<Targ>& x_left, const std::vector<Targ>& x_right);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_file_config(const std::string&                                   config_file_path,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                    const std::vector<Targ>& x_0, const Targ& delta);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_file_config(const std::string&                                   config_file_path,
                    const std::function<Tfit(const std::vector<Targ>&)>& fitness_function);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_string_config(const std::string&                                   config,
                      const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                      const std::vector<Targ>& x_left, const std::vector<Targ>& x_right,
                      const std::vector<Targ>& x_0);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_string_config(const std::string&                                   config,
                      const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                      const std::vector<Targ>& x_left, const std::vector<Targ>& x_right);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_string_config(const std::string&                                   config,
                      const std::function<Tfit(const std::vector<Targ>&)>& fitness_function,
                      const std::vector<Targ>& x_0, const Targ& delta);

    static std::shared_ptr<NOTKController<Targ, Tfit>>
    build_notk_string_config(const std::string&                                   config,
                      const std::function<Tfit(const std::vector<Targ>&)>& fitness_function);
};
}

#endif
