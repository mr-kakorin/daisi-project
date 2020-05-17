#ifndef NOTK_HEAVYBALLCONFIG_H
#define NOTK_HEAVYBALLCONFIG_H

#include "../../base/iconfig.h"

namespace notk
{
class HeavyBallConfig : public BaseOptConfig
{
  public:
    double              get_gradient_calc_step() const;
    double              get_weight() const;

    double gradient_calc_step;
    double weight;
};
}

SERIALIZIBLE_STRUCT(notk::HeavyBallConfig, srfl::CheckModes::FATAL,
                    SER_BASE()(double, gradient_calc_step, srfl::nan, 1e-16,
                               srfl::inf)(double, weight, srfl::nan, 0, srfl::inf))

#endif
