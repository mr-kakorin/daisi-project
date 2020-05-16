// #define _USE_MATH_DEFINES
// #define SEED 1
// #define BRNG VSL_BRNG_MCG31
#include "../base/AccelFlow.h"
#include "FlagStringsSolver.h"
#include "Tools.h"

#include <common_tools/constants.h>

#include <random>
// #include <mkl.h>

void AccelFlow::GenerateParticlesEnvelopes(double lambda, double z0, double Rchann)
{
    dynamicsTransv.clearandResize(nParticles);

    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    double P = gamma * beta;

    // VSLStreamStatePtr stream;

    // vslNewStream(&stream, BRNG, SEED);

    double R = channelRelRadius * Rchann;

    double emittanceXL = emittanceX * 1e-6;
    double emittanceYL = emittanceY * 1e-6;

    double dXmaxL = dXmax * beta;

    double betaX  = Xmax * Xmax / emittanceXL;
    double gammaX = dXmaxL * dXmaxL / emittanceXL;
    double betaY  = Xmax * Xmax / emittanceYL;
    double gammaY = dXmaxL * dXmaxL / emittanceYL;
    double AlphaX = sqrt(gammaX * betaX - 1);

    double phi = (180.0 / 3.14) * atan(2 * AlphaX / (gammaX - betaX)) / 2;

    if (dXmax < 0)
        AlphaX = -AlphaX;

    if (std::abs(gammaX * betaX - 1) < 1e-7)
        AlphaX = 0;

    double AlphaY = sqrt(gammaY * betaY - 1);
    if (dXmax < 0)
        AlphaY = -AlphaY;

    if (std::abs(gammaY * betaY - 1) < 1e-7)
        AlphaY = 0;

    for (int i = 0; i < nParticles; i++)
    {
        dynamicsTransv.S11x[i] = emittanceXL * betaX;
        dynamicsTransv.S12x[i] = -emittanceXL * AlphaX;
        dynamicsTransv.S22x[i] = emittanceXL * gammaX;

        dynamicsTransv.S11y[i]       = emittanceYL * betaY;
        dynamicsTransv.S12y[i]       = -emittanceYL * AlphaY;
        dynamicsTransv.S22y[i]       = emittanceYL * gammaY;
        dynamicsTransv.cellNumber[i] = -2;
    };

    // vdRngGaussian(0, stream, nParticles - 1, &dynamicsTransv.p[1], P, P * momentumSpread / 3.0);

    // if (phaseSpread == 0)
    // {
    //     vdRngUniform(0, stream, nParticles - 1, &dynamicsTransv.z[1], -lambda * beta + z0, z0);
    //     averagePhase = commtools::PI();
    // }
    // else
    //     vdRngGaussian(0, stream, nParticles - 1, &dynamicsTransv.z[1],
    //                   -averagePhase * lambda * beta / (2 * commtools::PI()) + z0,
    //                   phaseSpread / 3.0);

    // dynamicsTransv.p[i] = P;
    //	dynamicsTransv.z[i] = -lambda*beta + i*dZ + dZ / 2 + z0;

    dynamicsTransv.p[0] = P;
    dynamicsTransv.z[0] = -averagePhase * lambda * beta / (2 * commtools::PI()) + z0;

    dynamicsTransv.transvQ =
        impulseCurrent * lambda / (commtools::LIGHT_VELOCITY() * nParticlesTransv);

    if (nParticlesTransv == 1)
        dynamicsTransv.z[0] = z0;
};