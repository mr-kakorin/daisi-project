#ifndef NOTK_PARTICLESWARM_H
#define NOTK_PARTICLESWARM_H

#include <memory>

#include "../../base/ioptistep.h"

namespace notk
{
class ParticleSwarmConfig;

template <class Targ, class Tfit>
class ParticleSwarm final : public IOptimizationStep<Targ, Tfit>
{
  private:
    it_res_t<Targ, Tfit> do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                       const borders_t<Targ>& current_borders) override;

    it_res_t<Targ, Tfit> do_iteration(const size_t                iter_counter,
                                      const it_res_t<Targ, Tfit>& iter_result,
                                      const borders_t<Targ>& current_borders) override;

    borders_t<Targ> squeez_borders(const size_t iter_counter,
                                   const it_res_t<Targ, Tfit>& iter_result,
                                   const borders_t<Targ>& current_borders) override;

    std::shared_ptr<IOptimizationStep<Targ, Tfit>> m_searcher1d;

    //лучшие значения параметров для каждого агента алгоритма
    std::vector<std::vector<std::vector<Targ>>> bestAgentPopulations;
    //текушие значения параметров для каждого агента алгоритма
    std::vector<std::vector<std::vector<Targ>>> CurrentPopulations;
    //текушие скорости агентов
    std::vector<std::vector<std::vector<Targ>>> CurrentVelocities;
    //лучшие значения фитнесс-функции для каждого агента
    std::vector<std::vector<Tfit>> bestAgentsValues;
    //лучшее значение фитнесса для каждой нити
    std::vector<Tfit> bestThreadValues;
    //лучшее значение параметров для каждой нити
    std::vector<std::vector<Targ>> bestThreadPopulations;
    //лучшее параметров на итерации
    std::vector<Targ> bestIterationPopulations;
    Tfit              BestValue;
    volatile int      flagAllBreak;
    std::vector<Targ> bestPopulation;
    std::vector<Tfit> bestIterationValues;

    size_t nThreads;
    double phi1;
    double phi2;
    Tfit   lastBestValue;
    size_t iterRestart;
    size_t iterInsideRestart;
    size_t iterMaxRestart;
    bool   flagRestart;

  public:
    BUILD_CHILD(ParticleSwarm<Targ COMMA Tfit>, IOptimizationStep<Targ COMMA Tfit>,
                "Particle swarm")

    bool read_config(const boost::property_tree::ptree& config) override;
    std::shared_ptr<ParticleSwarmConfig> config() const;

    ParticleSwarm() = default;
};
}

#endif