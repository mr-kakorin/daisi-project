/*#define _USE_MATH_DEFINES
#define SEED 1
#include "LinacFlow.h"
#include <common_tools/constants.h>

// #define BRNG VSL_BRNG_MCG31
// 
// #include "mkl_vsl.h"

void LinacFlow::GenerateOneParticleRFQ(double lambda, double z0)
{
    oneParticle.resize(8);
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    double P       = gamma * beta;
    oneParticle[0] = P;
    oneParticle[1] = z0;

    double emittanceXL = emittanceX * 1e-6;
    double emittanceYL = emittanceY * 1e-6;

    double dXmaxL = dXmax * beta;

    double betaX  = Xmax * Xmax / emittanceXL;
    double gammaX = dXmaxL * dXmaxL / emittanceXL;
    double betaY  = Xmax * Xmax / emittanceYL;
    double gammaY = dXmaxL * dXmaxL / emittanceYL;
    double AlphaX = sqrt(gammaX * betaX - 1);
    if (dXmax < 0)
        AlphaX = -AlphaX;

    if (std::abs(gammaX * betaX - 1) < 1e-7)
        AlphaX = 0;

    double AlphaY = sqrt(gammaY * betaY - 1);
    if (dXmax < 0)
        AlphaY = -AlphaY;

    if (std::abs(gammaY * betaY - 1) < 1e-7)
        AlphaY = 0;

    oneParticle[2] = emittanceXL * betaX;
    oneParticle[3] = -emittanceXL * AlphaX;
    oneParticle[4] = emittanceXL * gammaX;

    oneParticle[5] = emittanceYL * betaY;
    oneParticle[6] = -emittanceYL * AlphaY;
    oneParticle[7] = emittanceYL * gammaY;

    oneParticleData.resize(4);
    oneParticleData[0] = -2;
}

std::vector<void*> LinacDynamicsLong::GetData()
{
    std::vector<void*> result(2);

    if (p.size() == 0)
        return result;

    result[0] = (void*)(&z[0]);
    result[1] = (void*)(&p[0]);

    return result;
}

double LinacFlow::GetOutVelocity() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(outEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta;
}

void LinacFlow::GenerateParticlesRFQ(double lambda, double z0, double Rchann)
{
    dynamicsLong.clearandResize(nParticlesLong);
    dynamicsTransv.clearandResize(nParticlesTransv);

    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    double P = gamma * beta;

    // VSLStreamStatePtr stream;

    // vslNewStream(&stream, BRNG, SEED);

    // vdRngGaussian(0, stream, nParticlesLong - 1, &dynamicsLong.p[1], P, P * momentumSpread / 3);

    // vdRngUniform(0, stream, nParticlesLong - 1, &dynamicsLong.z[1], -lambda * beta + z0, z0);

    dynamicsLong.p[0] = P;
    dynamicsLong.z[0] = -lambda * beta / 2 + z0;

    double R           = channelRelRadius * Rchann;
    dynamicsLong.longQ = impulseCurrent * lambda /
                         (commtools::LIGHT_VELOCITY() * nParticlesLong * commtools::PI() * R * R);

    if (nParticlesTransv == 1)
        dynamicsLong.z[0] = z0;

    for (int i = 0; i < nParticlesLong; i++)
    {
        dynamicsLong.cellNumber[i] = -2;
    };

    double emittanceXL = emittanceX * 1e-6;
    double emittanceYL = emittanceY * 1e-6;

    double dXmaxL = dXmax * beta;

    double betaX  = Xmax * Xmax / emittanceXL;
    double gammaX = dXmaxL * dXmaxL / emittanceXL;
    double betaY  = Xmax * Xmax / emittanceYL;
    double gammaY = dXmaxL * dXmaxL / emittanceYL;
    double AlphaX = sqrt(gammaX * betaX - 1);
    if (dXmax < 0)
        AlphaX = -AlphaX;

    if (std::abs(gammaX * betaX - 1) < 1e-7)
        AlphaX = 0;

    double AlphaY = sqrt(gammaY * betaY - 1);
    if (dXmax < 0)
        AlphaY = -AlphaY;

    if (std::abs(gammaY * betaY - 1) < 1e-7)
        AlphaY = 0;

    double dZ = lambda * beta / nParticlesTransv;

    for (int i = 0; i < nParticlesTransv; i++)
    {
        dynamicsTransv.S11x[i] = emittanceXL * betaX;
        dynamicsTransv.S12x[i] = -emittanceXL * AlphaX;
        dynamicsTransv.S22x[i] = emittanceXL * gammaX;

        dynamicsTransv.S11y[i]       = emittanceYL * betaY;
        dynamicsTransv.S12y[i]       = -emittanceYL * AlphaY;
        dynamicsTransv.S22y[i]       = emittanceYL * gammaY;
        dynamicsTransv.cellNumber[i] = -2;
        dynamicsTransv.z[i]          = -lambda * beta + i * dZ + dZ / 2 + z0;
    };

    // vdRngGaussian(0, stream, nParticlesTransv - 1, &dynamicsTransv.p[1], P, P * momentumSpread /
    // 3);

    // vdRngUniform(0, stream, nParticlesTransv - 1, &dynamicsTransv.z[1], -lambda*beta + z0, z0);
    // dynamicsTransv.p[i] = P;
    //	dynamicsTransv.z[i] = -lambda*beta + i*dZ + dZ / 2 + z0;

    dynamicsTransv.p[0] = P;
    dynamicsTransv.z[0] = -lambda * beta / 2 + z0;

    dynamicsTransv.transvQ =
        impulseCurrent * lambda / (commtools::LIGHT_VELOCITY() * nParticlesTransv);

    if (nParticlesTransv == 1)
        dynamicsTransv.z[0] = z0;
};
double LinacFlow::GetStartVelocity() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta;
};
double LinacFlow::GetStartMomentum() const
{
    double restEnergy = -mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() /
                        commtools::ELECTRON_CHARGE(); // ������� ����� ������� � ���������������
    double gamma = (restEnergy + std::abs(averageEnergy)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta * gamma;
}
double LinacFlow::GetAlpha() const
{
    return charge / (mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY());
};
double LinacFlow::GetRestEnergy() const
{
    return mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY();
};

void LinacDynamicsLong::clearandResize(int nParticles)
{
    NparticlesLong = nParticles;
    p.clear();
    z.clear();

    Ez.clear();
    cellNumber.clear();

    p.resize(nParticles);
    z.resize(nParticles);
    Ez.resize(nParticles);
    cellNumber.resize(nParticles);
};

void LinacDynamicsTransv::clearandResize(int nParticles)
{
    NparticlesTransv = nParticles;
    p.clear();
    z.clear();
    S11x.clear();
    S12x.clear();
    S22x.clear();
    S11y.clear();
    S12y.clear();
    S22y.clear();
    Ex.clear();
    Ey.clear();
    Ez.clear();
    cellNumber.clear();

    p.resize(nParticles);
    z.resize(nParticles);
    S11x.resize(nParticles);
    S12x.resize(nParticles);
    S22x.resize(nParticles);
    S11y.resize(nParticles);
    S12y.resize(nParticles);
    S22y.resize(nParticles);
    Ex.resize(nParticles);
    Ey.resize(nParticles);
    Ez.resize(nParticles);
    EzCol.resize(nParticles);
    cellNumber.resize(nParticles);
};

int LinacDynamicsTransv::checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                                           const std::vector<double>& LL,
                                           const std::vector<double>& MinumalRadii, double lambda,
                                           double relChann)
{
    double zAv = 0;
    double pAv = 0;

    int    n    = 0;
    double maxR = 0;
    double Lend = LL.back();
    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    for (int i = 0; i < z.size(); i++)
    {
        if (cellNumber[i] != -1)
            n++;
    };
    for (int i = 0; i < z.size(); i++)
    {
        //	if (cellNumber[i] != -1)
        //	{
        double R = std::max(S11x[i], S11y[i]);
        if (R > maxR)
            maxR = R;
        zAv      = zAv + z[i];
        pAv      = pAv + p[i];
        //	};
    };
    zAv = zAv / z.size();
    pAv = pAv / z.size();

    cell = 0;
    for (cell = 0; cell < LL.size(); cell++)
    {
        if (LL[cell] > z[0])
        {
            cell--;
            break;
        }
    };
    if (cell == -1)
        cell = 0;
    if (cell >= LL.size() - 1)
        cell = LL.size() - 2;

    double Rchannel = MinumalRadii[cell];
    OutParams[0].push_back(Rchannel);
    OutParams[1].push_back(sqrt(maxR));

    double transmission = 0;
    double acceleration = 0;

    for (int i = 0; i < z.size(); i++)
    {
        // if (z[i]>Lend || z[i]<0)
        //	cellNumber[i] = -1;

        if (z[i] > Lend)
        {
            if (cellNumber[i] != -1)
            {
                OutParams[8].push_back(S11x[i]);
                OutParams[9].push_back(S12x[i]);
                OutParams[10].push_back(S22x[i]);
                OutParams[11].push_back(S11y[i]);
                OutParams[12].push_back(S12y[i]);
                OutParams[13].push_back(S22y[i]);
            }
            cellNumber[i] = -1;
        }
        double R = sqrt(S11x[i]);
        if (S11x[i] < 1e-8)
            R = 1.0000e-04;

        if (R < Rchannel * relChann || cell < 20)
            transmission++;
        else
            cellNumber[i] = -1;

        double dPh = 2 * commtools::PI() * (zAv - z[i]) / (pAv * lambda);

        if (dPh < commtools::PI())
            acceleration++;

        //	if (std::abs(pAv - p[i]) / pAv < 0.04)
        //	acceleration++;
    }

    if (acceleration == 0)
    {
        int tt = 0;
    }
    OutParams[3].push_back(acceleration / z.size());

    OutParams[2].push_back(transmission / z.size());
    OutParams[4].push_back(pAv);

    return n;
};

int LinacDynamicsLong::checkBadParticles(int& cell, std::vector<std::vector<float>>& OutParams,
                                         const std::vector<double>& LL,
                                         const std::vector<double>& MinumalRadii, double lambda)
{
    double zAv = z[0];
    double pAv = 0;

    int    n    = 0;
    double maxR = 0;
    double Lend = LL.back();
    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    for (int i = 0; i < z.size(); i++)
    {
        if (cellNumber[i] != -1)
            n++;
    };
    for (int i = 0; i < z.size(); i++)
    {
        //		zAv = zAv + z[i];
        pAv = pAv + p[i];
    };
    //	zAv = zAv / z.size();
    pAv = pAv / z.size();

    cell = 0;
    for (cell = 0; cell < LL.size(); cell++)
    {
        if (LL[cell] > zAv)
        {
            cell--;
            break;
        }
    };
    if (cell == -1)
        cell = 0;
    if (cell >= LL.size() - 1)
        cell = LL.size() - 2;

    double transmission = 0;
    double acceleration = 0;

    for (int i = 0; i < z.size(); i++)
    {
        //	if (z[i]>Lend || z[i]<0)
        //		cellNumber[i] = -1;

        if (z[i] > Lend)
            cellNumber[i] = -1;
        double dPh        = 2 * commtools::PI() * (zAv - z[i]) / (pAv * lambda);

        if (dPh < commtools::PI())
            acceleration++;

        //	if (std::abs(pAv - p[i]) / pAv < 0.04)
        //	acceleration++;
    }

    if (acceleration == 0)
    {
        int tt = 0;
    }

    OutParams[3].push_back(acceleration / z.size());

    OutParams[4].push_back(pAv);

    return n;
};
*/