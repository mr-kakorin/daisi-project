#include "heavyballconfig.h"

namespace notk
{
double HeavyBallConfig::get_gradient_calc_step() const
{
    return gradient_calc_step;
}
double HeavyBallConfig::get_weight() const
{
    return weight;
}
}