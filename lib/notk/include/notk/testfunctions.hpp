#ifndef NOTK_TESTFUNCTIONS_H
#define NOTK_TESTFUNCTIONS_H

#include <vector>

namespace notk
{
double Rosenbrock(const std::vector<double>& population);
double Rastrigin(const std::vector<double>& population);
double Sphere(const std::vector<double>& population);
double Ackley(const std::vector<double>& population);
}

#endif
