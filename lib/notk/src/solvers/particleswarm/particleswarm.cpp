#include <numeric>
#include <random>

#include <serreflection/read_json.hpp>

#include "particleswarm.h"
#include "particleswarmconfig.h"

namespace notk
{

template class ParticleSwarm<double, double>;
template class ParticleSwarm<float, float>;

REGISTER_CHILD(ParticleSwarm<double COMMA double>, IOptimizationStep<double COMMA double>,
               "Particle swarm")
REGISTER_CHILD(ParticleSwarm<float COMMA float>, IOptimizationStep<float COMMA float>,
               "Particle swarm")

template <class Targ, class Tfit>
it_res_t<Targ, Tfit>
ParticleSwarm<Targ, Tfit>::do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                         const borders_t<Targ>& current_borders)
{
    nThreads           = 1;
    size_t problem_dim = current_borders.first.size();

    bestAgentPopulations.resize(nThreads);
    CurrentPopulations.resize(nThreads);
    bestAgentsValues.resize(nThreads);
    CurrentVelocities.resize(nThreads);

    size_t nAgentsPErThread = config()->get_n_agents() / nThreads + 1;

    for (size_t thread = 0; thread < nThreads; thread++)
    {
        bestAgentPopulations[thread].resize(nAgentsPErThread);
        CurrentPopulations[thread].resize(nAgentsPErThread - 1);
        bestAgentsValues[thread].resize(nAgentsPErThread);
        CurrentVelocities[thread].resize(nAgentsPErThread);
    }

    //выбираем весовые параметры в методе роя случайно
    phi1              = double(rand()) / RAND_MAX;
    phi2              = 1 - phi1;
    iterRestart       = 0;
    iterInsideRestart = 0;
    iterMaxRestart    = 100;
    flagRestart       = false;

    //генерируем начальную популяцию

    std::default_random_engine           generator;
    std::normal_distribution<Targ>       distribution(0, 1);
    std::uniform_real_distribution<Targ> distributionUn(0, 1);

    for (size_t thread = 0; thread < nThreads; thread++)
    {
        for (size_t j = 0; j < CurrentPopulations[thread].size(); j++)
        {
            CurrentPopulations[thread][j].resize(problem_dim);

            for (size_t k = 0; k < problem_dim; k++)
            {
                if (std::abs(current_borders.second[k] - current_borders.first[k]) < 1e-15)
                {
                    CurrentPopulations[thread][j][k] = 0;
                }
                else
                {
                    if (config()->get_is_gauss())
                    {
                        CurrentPopulations[thread][j][k] =
                            current_borders.first[k] +
                            distributionUn(generator) * (2 * current_borders.second[k]);
                    }
                    else
                    {
                        auto sigma = (current_borders.second[k] - current_borders.first[k]) / 6;
                        auto mu    = (current_borders.second[k] + current_borders.first[k]) / 2;
                        CurrentPopulations[thread][j][k] = sigma * distribution(generator) + mu;
                    }
                }
            }

            bestAgentPopulations[thread][j] = CurrentPopulations[thread][j];
            bestAgentsValues[thread][j]     = std::numeric_limits<Tfit>::max();
            CurrentVelocities[thread][j].resize(CurrentPopulations[thread][j].size());
        }
        CurrentPopulations[thread].push_back(iter_result.first);
        CurrentVelocities[thread].back().resize(CurrentPopulations[thread].back().size());
        bestAgentPopulations[thread].back() = CurrentPopulations[thread].back();
    }

    //лучшее значение фитнесса для каждой нити
    bestThreadValues.resize(nThreads);
    //лучшее значение параметров для каждой нити
    bestThreadPopulations.resize(nThreads);
    //лучшее параметров на итерации
    BestValue    = 1e20;
    flagAllBreak = 0;

    return iter_result;
}

template <class Targ, class Tfit>
it_res_t<Targ, Tfit>
ParticleSwarm<Targ, Tfit>::do_iteration(const size_t                iter_counter,
                                        const it_res_t<Targ, Tfit>& iter_result,
                                        const borders_t<Targ>& current_borders)
{
    const size_t problem_dim = current_borders.first.size();

    volatile size_t thread   = 0;
    bestThreadValues[thread] = std::numeric_limits<Tfit>::max();

    for (size_t j = 0; j < CurrentPopulations[thread].size(); j++)
    {
        // в параллельном режиме вычислем фитнесс для всех агентов
        double f = this->fitness(CurrentPopulations[thread][j]);
        // лучшая точка для каждого отдельного агента
        if (f < bestAgentsValues[thread][j])
        {
            bestAgentsValues[thread][j]     = f;
            bestAgentPopulations[thread][j] = CurrentPopulations[thread][j];
        }
        // лучшая точка для каждой нити
        if (f < bestThreadValues[thread])
        {
            bestThreadValues[thread]      = f;
            bestThreadPopulations[thread] = CurrentPopulations[thread][j];
        }
    }

    //лучшее значение на итерации
    double bestIterationValue = 1e20;
    for (size_t thread = 0; thread < nThreads; thread++)
    {
        if (bestThreadValues[thread] < bestIterationValue)
        {
            bestIterationValue       = bestThreadValues[thread];
            bestIterationPopulations = bestThreadPopulations[thread];
        }
    }

    bestIterationValues.push_back(bestIterationValue);

    //лучшее значение вообще
    if (bestIterationValue < BestValue)
    {
        BestValue      = bestIterationValue;
        bestPopulation = bestIterationPopulations;
    }

    //вывод

    std::stringstream ss2;
    ss2 << std::scientific << BestValue;

    it_res_t<Targ, Tfit> result = std::make_pair(bestPopulation, this->fitness(bestPopulation));

    //обновление скоростей и положений агентов
    for (size_t thread = 0; thread < nThreads; thread++)
    {
        for (size_t j = 0; j < CurrentPopulations[thread].size(); j++)
        {
            int nn1 = rand() % 100;
            int nn2 = rand() % 100;

            double r1 = double(nn1) / 100;
            double r2 = double(nn2) / 100;

            for (size_t i = 0; i < CurrentPopulations[thread][j].size(); i++)
            {
                CurrentVelocities[thread][j][i] =
                    CurrentVelocities[thread][j][i] +
                    phi1 * r1 *
                        (bestAgentPopulations[thread][j][i] - CurrentPopulations[thread][j][i]) +
                    +phi2 * r2 * (bestPopulation[i] - CurrentPopulations[thread][j][i]);

                CurrentPopulations[thread][j][i] =
                    CurrentPopulations[thread][j][i] + CurrentVelocities[thread][j][i];
            }
            for (size_t k = 0; k < problem_dim; k++)
            {
                if (CurrentPopulations[thread][j][k] > current_borders.second[k])
                    CurrentPopulations[thread][j][k] = current_borders.second[k];

                if (CurrentPopulations[thread][j][k] < current_borders.first[k])
                    CurrentPopulations[thread][j][k] = current_borders.first[k];
            }
        }
    }

    // если алгоритм стопорится на одном лучшем значении более чем iterMaxRestart итераций,
    // то переразыгрываем весовые значения phi1 и phi2 если выполняется условие как в алгоритме
    // имитации
    // отжига с вероятностью p = exp(-double(iterRestart - iterMaxRestart) /
    // double(iterMaxRestart));
    if (iter_counter > 1 && lastBestValue > BestValue)
        flagRestart = false;

    if (iter_counter > 1 && lastBestValue == BestValue)
        iterRestart++;
    else
        iterRestart = 0;

    lastBestValue = BestValue;

    double p = 1;
    if (iterRestart > iterMaxRestart)
    {
        p = exp(-double(iterRestart - iterMaxRestart) / double(iterMaxRestart));
    }
    double r = double(rand()) / RAND_MAX;
    if (r > p)
    {

        phi1        = double(rand()) / RAND_MAX;
        phi2        = 1 - phi1;
        iterRestart = 0;
    }

    return result;
}

template <class Targ, class Tfit>
borders_t<Targ> ParticleSwarm<Targ, Tfit>::squeez_borders(const size_t                iter_counter,
                                                          const it_res_t<Targ, Tfit>& iter_result,
                                                          const borders_t<Targ>& current_borders)
{
    return current_borders;
}
template <class Targ, class Tfit>
bool ParticleSwarm<Targ, Tfit>::read_config(const boost::property_tree::ptree& config)
{
    this->m_config = srfl::read_json<ParticleSwarmConfig>(config);

    return this->m_config ? true : false;
}
template <class Targ, class Tfit>
std::shared_ptr<ParticleSwarmConfig> ParticleSwarm<Targ, Tfit>::config() const
{
    return std::static_pointer_cast<ParticleSwarmConfig>(this->m_config);
}
}