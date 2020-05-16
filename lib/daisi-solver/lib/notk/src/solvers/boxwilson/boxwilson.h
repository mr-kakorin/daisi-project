#ifndef NOTK_BOXWILSON_H
#define NOTK_BOXWILSON_H

#include <memory>

#include "../../base/ioptistep.h"
#include "boxwilsonconfig.h"

namespace notk
{

template <class Targ, class Tfit>
class BoxWilson final : public IOptimizationStep<Targ, Tfit>
{
  private:
    it_res_t<Targ, Tfit> do_iteration(const size_t                iter_counter,
                                      const it_res_t<Targ, Tfit>& iter_result,
                                      const borders_t<Targ>& current_borders) override;

    borders_t<Targ> squeez_borders(const size_t iter_counter,
                                   const it_res_t<Targ, Tfit>& iter_result,
                                   const borders_t<Targ>& current_borders) override;

  public:
    BUILD_CHILD(BoxWilson<Targ COMMA Tfit>, IOptimizationStep<Targ COMMA Tfit>, "Box-Wilson")
    BoxWilson() = default;
    bool read_config(const boost::property_tree::ptree& config) override;

    std::shared_ptr<BoxWilsonConfig> config() const;
};
}

#endif