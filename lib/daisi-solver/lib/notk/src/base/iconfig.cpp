#include <cmath>

#include <boostlog0/boostlog.h>
#include <common_tools/arrays.h>

#include "iconfig.h"

namespace notk
{

size_t BaseOptConfig::get_maxNumIterations() const
{
    return maximal_iterations;
}

template bool
BaseOptConfig::check_accuracy<double, double>(const it_res_t<double, double>& result_prev,
                                              const it_res_t<double, double>& result_cur,
                                              const bool do_log, const bool is_first) const;
template bool BaseOptConfig::check_accuracy<float, float>(const it_res_t<float, float>& result_prev,
                                                          const it_res_t<float, float>& result_cur,
                                                          const bool do_log,
                                                          const bool is_first) const;

template <class Targ, class Tfit>
bool BaseOptConfig::check_accuracy(const it_res_t<Targ, Tfit>& result_prev,
                                   const it_res_t<Targ, Tfit>& result_cur, const bool do_log,
                                   const bool is_first) const
{
    double currentAccuracy;
    auto calcAccuracy = [&](const std::vector<Targ>& x1, const std::vector<Targ>& x2) {
        // вычисление точности
        //     std::vector<Targ>
        auto difference = commtools::apply_op(x2, x1, commtools::Operations::MINUS);

        switch (accuracy_type)
        {
        // по абсолютным значениям
        case AccuracyType::ABSOLUTE:
            std::for_each(difference.begin(), difference.end(),
                          [](auto& val) { val = std::abs(val); });
            break;
        // по относительным значениям
        case AccuracyType::RELATIVE:
            for (size_t i = 0; i < difference.size(); i++)
            {
                if (std::abs(x2[i]) < 1e-9)
                {
                    difference[i] = std::abs(difference[i]);
                }
                else
                {
                    difference[i] = std::abs(difference[i] / x2[i]);
                }
            }
            break;
        case AccuracyType::ZERO:
            difference = x2;
            std::for_each(difference.begin(), difference.end(),
                          [](auto& val) { val = std::abs(val); });
            break;
        default:
            throw std::runtime_error("getAccuracy::unexpected accuracyType");
        }

        auto result = std::max_element(difference.begin(), difference.end());
        return *result;
    };

    switch (accuracy_source_type)
    {
    // вычисление точности по аргументу фитнеса
    case AccuracySourceType::ARGUMENT:
        currentAccuracy = calcAccuracy(result_prev.first, result_cur.first);
        break;
    // вычисление точности по фитнесу
    case AccuracySourceType::FITNESS:
        currentAccuracy = calcAccuracy(std::vector<Targ>{result_prev.second},
                                       std::vector<Targ>{result_cur.second});
        break;
    default:
        throw std::runtime_error("getAccuracy::unexpected accuracySourceType");
    }
    if (do_log)
    {
        LOG(sev_lvl::trace) << "call: check_accuracy done, currentAccuracy = " << currentAccuracy;
    }
    if (accuracy_source_type == AccuracySourceType::FITNESS && accuracy_type == AccuracyType::ZERO)
    {
        return currentAccuracy < required_accuracy;
    }
    return currentAccuracy < required_accuracy && !is_first;
}
}
