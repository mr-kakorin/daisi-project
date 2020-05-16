#ifndef NOTK_HEAVYBALL_H
#define NOTK_HEAVYBALL_H

#include <memory>

#include "../../base/ioptistep.h"

namespace notk
{
class HeavyBallConfig;

template <class Targ, class Tfit>
class HeavyBall : public IOptimizationStep<Targ, Tfit>
{
  private:
    virtual it_res_t<Targ, Tfit>
    do_iteration(const size_t iter_counter, const it_res_t<Targ, Tfit>& iter_result,
                 const borders_t<Targ>& current_borders) override final;

    virtual borders_t<Targ>
    squeez_borders(const size_t iter_counter, const it_res_t<Targ, Tfit>& iter_result,
                   const borders_t<Targ>& current_borders) override final;

    std::vector<Targ> last_x;

    std::shared_ptr<IOptimizationStep<Targ, Tfit>> m_searcher1d;

  public:
    BUILD_CHILD(HeavyBall<Targ COMMA Tfit>, IOptimizationStep<Targ COMMA Tfit>, "Heavy ball")

    HeavyBall() = default;
    virtual bool read_config(const boost::property_tree::ptree& config) override final;
    std::shared_ptr<HeavyBallConfig> config() const;
};
}

#endif