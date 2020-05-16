#include "EmitterDevice2daxsVectorized.h"
#include <common_tools/constants.h>
#include <cmath>

#include "mkl_vml.h"

#include <vector>
std::vector<float> steps(512);
std::vector<float> tmp(512);
std::vector<float> tmp1(512);

std::vector<std::vector<float>> stepsT(10);
std::vector<std::vector<float>> tmpT(10);
std::vector<std::vector<float>> tmp1T(10);

double IntegrateEnergy(double a, double b, double t0)
{
    int    s    = 300;
    double step = (b - a) / s;

    double x = a;
    steps[0] = a;

    for (int i = 1; i < s + 1; i++)
    {
        steps[i] = steps[i - 1] + step;
    };

    for (int i = 0; i < s + 1; i++)
    {
        tmp[i] = (2 / t0) * sqrt(steps[i] / (commtools::PI() * t0)) * exp(-steps[i] / t0);
    }

    double result = 0;

    for (int i = 0; i < s; i++)
    {
        result = result + tmp[i] + tmp[i + 1];
    };
    result = result * step / 2;

    return result;
};

double IntegrateVelocity(double a, double b, double t0, double m, double charge, int thread)
{
    int    s    = 300;
    double step = (b - a) / s;

    double x = a;
    stepsT[thread].resize(512);
    tmpT[thread].resize(512);
    stepsT[thread][0] = a;

    for (int i = 1; i < s + 1; i++)
    {
        stepsT[thread][i] = stepsT[thread][i - 1] + step;
    };

    for (int i = 0; i < s + 1; i++)
    {

        tmpT[thread][i] = sqrt(m / (2 * commtools::PI() * t0 * std::abs(charge))) *
                          exp(-(m * stepsT[thread][i] * stepsT[thread][i]) / (2 * t0 * std::abs(charge)));
    }

    double result = 0;

    for (int i = 0; i < s; i++)
    {
        result = result + tmpT[thread][i] + tmpT[thread][i + 1];
    };
    result = result * step / 2;

    return result;
};

float IntegrateCurrent(float a, float b, float t0)
{
    int   s      = int((b - a) * 200 / commtools::PI());
    float result = 0;
    float step   = (b - a) / s;
    // double f1 = f(a, param);
    //	double f1 = exp(-(std::sin(a)*std::sin(a)) / (std::sin(t0)*std::sin(t0)));
    //	double currentP = a;
    //	double f2;

    steps[0] = a;

    for (int i = 1; i < s + 1; i++)
    {
        steps[i] = steps[i - 1] + step;
    };

    float st = std::sin(t0) * std::sin(t0);
    vsSin(s + 1, &steps[0], &tmp[0]);

#pragma ivdep
    for (int i = 0; i < s + 1; i++)
    {
        tmp[i] = -tmp[i] * tmp[i] / (st);
    };

    vsExp(s + 1, &tmp[0], &tmp1[0]);

    for (int i = 0; i < s; i++)
    {
        result = result + tmp1[i] + tmp1[i + 1];
    };
    result = result * step / 2;
    if (std::isinf(result))
    {
        int tt = 0;
    };
    return result;
};

template void EmitterDevice2daxsVectorized<float>(float* a, float* b);
template void EmitterDevice2daxsVectorized<double>(double* a, double* b);

template <class PointType>
void EmitterDevice2daxsVectorized(PointType* rPointer, PointType* zPointer, PointType* pphiPointer,
                                  PointType* prPointer, PointType* phiPointer, PointType* pzPointer,
                                  PointType* qPointer, PointType* cellsNumbersPointer, PointType dt)
{
    PointType dphPhiR;
    PointType phRPhi1;
    PointType phPhiR;
    PointType current;
    PointType phRPhi0;
    PointType index;
#pragma ivdep
#pragma vector always
    for (int i3 = 0; i3 < nParticlesPhiR; i3++)
    {
        /*dphPhiR = DPH1 + i3*KDPH;
        phRPhi1 = phRPhi0 + dphPhiR;
        phPhiR = (phRPhi0 + phRPhi1) / 2;

        if (nParticlesPhiR == 1)
        phPhiR = 0;

        current = 2 * CphRPhi * Dmath::integral(phRPhi0, phRPhi1, andgleDistribution, phiPhiR);

        current = 1;
        phRPhi0 = phRPhi1;
        if (k < empty)
        index = EmptyPlaces[k];
        else
        {
        index = nowParticles + k1;
        k1++;
        }*/

        dphPhiR = 0;
        phRPhi1 = 0;
        phPhiR  = 0;

        phPhiR = 0;

        current = 1;

        current = 1;
        phRPhi0 = phRPhi1;

        index = i3;

        rPointer[index] = r;
        zPointer[index] = z;

        PointType prz              = pTotal * std::cos(phPhiR);
        pphiPointer[index]         = pTotal * std::sin(phPhiR) * r;
        prPointer[index]           = -prz * std::cos(phRZ) * sign;
        phiPointer[index]          = 0;
        pzPointer[index]           = prz * std::sin(phRZ);
        qPointer[index]            = current * dt * chargeSign;
        cellsNumbersPointer[index] = 0;
        //	current = particlesData->q[nowParticles + k];
        k++;
    };
};
