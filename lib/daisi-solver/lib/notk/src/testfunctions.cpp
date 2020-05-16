#include <cmath>
#include <cstddef>

#include <notk/testfunctions.hpp>

namespace notk
{
double Rosenbrock(const std::vector<double>& population)
{
    double f = 0;
    for (size_t i = 0; i < population.size() - 1; i++)
    {
        f = f + (100.0 * (population[i + 1] - population[i] * population[i]) *
                     (population[i + 1] - population[i] * population[i]) +
                 (population[i] - 1) * (population[i] - 1));
    }
    return f;
}

double Rastrigin(const std::vector<double>& population)
{
    double f = 0;
    double A = 10;
    for (size_t i = 0; i < population.size(); i++)
    {
        f = f + (population[i] * population[i] - A * std::cos(2 * 3.14 * population[i]));
    }
    f = A * population.size() + f;

    return f;
}
double Sphere(const std::vector<double>& population)
{
    double f = 0;

    for (size_t i = 0; i < population.size(); i++)
    {
        f = f + (population[i] * population[i]);
    }

    return f;
}
double Ackley(const std::vector<double>& population)
{
    double f  = 0;
    double f1 = 0;

    double a = 20;
    double b = 0.2;
    double c = 2 * 3.14159265359;

    for (size_t i = 0; i < population.size(); i++)
    {
        f = f + (population[i] * population[i]);
    }
    f = -b * std::sqrt(f / static_cast<double>(population.size()));

    for (size_t i = 0; i < population.size(); i++)
    {
        f1 = f1 + std::cos(c * population[i]);
    }
    f1 = f1 / static_cast<double>(population.size());

    return -a * std::exp(f) - std::exp(f1) + 20.0 + std::exp(1.0);
}
}