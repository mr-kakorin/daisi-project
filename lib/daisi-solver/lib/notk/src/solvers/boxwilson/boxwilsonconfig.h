#ifndef NOTK_BOXWILSON_CONFIG_H
#define NOTK_BOXWILSON_CONFIG_H

#include "../../base/iconfig.h"

namespace notk
{
class BoxWilsonConfig final : public BaseOptConfig
{
  public:
    double get_gradient_calc_step() const;

  private:
    double gradient_calc_step;
};
}

#endif
