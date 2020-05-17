#ifndef NOTK_GAUSSEARCH_H
#define NOTK_GAUSSEARCH_H

#include <memory>

#include "../../base/ioptistep.h"

namespace notk
{
class GaussSearchConfig;

template <class Targ, class Tfit>
class GaussSearch final : public IOptimizationStep<Targ, Tfit>
{
  private:
    Tfit fitness1dWrapper(const double x, const std::vector<Targ>& all, const size_t dim) const;

    Tfit single_iteration(const std::vector<size_t>& indexes, std::vector<Targ>& x_current,
                          const borders_t<Targ>& current_borders) const;

    it_res_t<Targ, Tfit> do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                       const borders_t<Targ>& current_borders) override;

    it_res_t<Targ, Tfit> do_iteration(const size_t                iter_counter,
                                      const it_res_t<Targ, Tfit>& iter_result,
                                      const borders_t<Targ>& current_borders) override;

    borders_t<Targ> squeez_borders(const size_t iter_counter,
                                   const it_res_t<Targ, Tfit>& iter_result,
                                   const borders_t<Targ>& current_borders) override;

    std::shared_ptr<IOptimizationStep<Targ, Tfit>> m_searcher1d;

  public:
    BUILD_CHILD(GaussSearch<Targ COMMA Tfit>, IOptimizationStep<Targ COMMA Tfit>)

    GaussSearch() = default;
    bool read_config(const boost::property_tree::ptree& config) override;
    std::shared_ptr<GaussSearchConfig> config() const;
};
}

#endif
