#ifndef NOTK_I_OPTIMCONFIG_H
#define NOTK_I_OPTIMCONFIG_H

#include <serreflection/defines.hpp>

#include <notk/types.hpp>

#include "../enums.h"

namespace notk
{
class BaseOptConfig
{
  public:
    size_t get_maxNumIterations() const;

    template <class Targ, class Tfit>
    bool check_accuracy(const it_res_t<Targ, Tfit>& result_prev,
                        const it_res_t<Targ, Tfit>& result_cur, const bool do_log,
                        const bool is_first) const;

    int maximal_iterations;

    AccuracyType       accuracy_type;
    AccuracySourceType accuracy_source_type;
    double             required_accuracy;
};
}

#define SER_BASE()                                                                                 \
    (int, maximal_iterations, srfl::nan, 1,                                                        \
     srfl::inf)(notk::AccuracyType, accuracy_type,                                                 \
                DEF_D())(notk::AccuracySourceType, accuracy_source_type,                           \
                         DEF_D())(double, required_accuracy, srfl::nan, 0, srfl::inf)

#endif
