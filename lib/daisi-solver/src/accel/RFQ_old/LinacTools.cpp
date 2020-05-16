#include <random>

#ifdef WIN32
#include <common_tools/include/constants.h>
#endif
#ifndef WIN32
#include <common_tools/constants.h>
#endif

#include "../base/LinacTools.h"
#include "Results.h"

std::default_random_engine       generator_lt;
std::normal_distribution<double> distribution_lt(0, 1);
// std::vector<std::vector<std::vector<double>>> F;
std::vector<std::vector<double>> Fn;
std::vector<std::vector<double>> FnOnePartilce;
std::vector<double>              rightV;
std::vector<double>              Vtmp;

void solveMatrix(int n, std::vector<double>& a, std::vector<double> c, std::vector<double>& b,
                 std::vector<double> f, std::vector<double>& x)
{
    double m;
    for (int i = 1; i < n; i++)
    {
        m    = a[i - 1] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    x[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--)
        x[i]   = (f[i] - b[i] * x[i + 1]) / c[i];
}

void LinacTools::chargeCalculate(LinacDynamicsTransv& dyn, std::vector<double>& W1,
                                 std::vector<double>& rho, double r)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        int cell = dyn.cellNumber[i];
        if (cell == -1)
            continue;

        double rx = sqrt(dyn.S11x[i]);
        double ry = sqrt(dyn.S11y[i]);

        double l = 1e-8;

        if (dyn.S11x[i] < l)
            rx = sqrt(l);

        if (dyn.S11y[i] < l)
            ry = sqrt(l);

        double rr     = r * r;
        rho[cell]     = rho[cell] + dyn.transvQ * W1[i] / (commtools::PI() * rr);
        rho[cell + 1] = rho[cell + 1] + dyn.transvQ * (1 - W1[i]) / (commtools::PI() * rr);
    }
};

void LinacTools::calculateF(std::vector<double>& dyn, std::vector<double>& oneParticleData,
                            double alpha, std::vector<double>& f)
{

    /*oneParticle[2] = emittanceXL*betaX;
    oneParticle[3] = -emittanceXL*AlphaX;
    oneParticle[4] = emittanceXL*gammaX;

    oneParticle[5] = emittanceYL*betaY;
    oneParticle[6] = -emittanceYL*AlphaY;
    oneParticle[7] = emittanceYL*gammaY;*/

    f[0] = alpha * oneParticleData[1];
    f[1] = dyn[0];

    f[2] = 2 * dyn[3];
    f[3] = oneParticleData[2] * alpha * dyn[2] + dyn[4];
    f[4] = 2 * oneParticleData[2] * alpha * dyn[3];

    f[5] = 2 * dyn[6];
    f[6] = oneParticleData[3] * alpha * dyn[5] + dyn[7];
    f[7] = 2 * oneParticleData[3] * alpha * dyn[6];
}

void LinacTools::calculateF(LinacDynamicsTransv& dyn, double alpha,
                            std::vector<std::vector<double>>& f)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        f[0][i] = alpha * dyn.Ez[i];
        f[1][i] = dyn.p[i];
        f[2][i] = 2 * dyn.S12x[i];
        f[3][i] = dyn.Ex[i] * alpha * dyn.S11x[i] + dyn.S22x[i];
        f[4][i] = 2 * dyn.Ex[i] * alpha * dyn.S12x[i];
        f[5][i] = 2 * dyn.S12y[i];
        f[6][i] = dyn.Ey[i] * alpha * dyn.S11y[i] + dyn.S22y[i];
        f[7][i] = 2 * dyn.Ey[i] * alpha * dyn.S12y[i];
    }
}
void LinacTools::updateMomentumsAndPositions(std::vector<LinacDynamicsTransv>& dyn, double dt,
                                             double alpha, int step,
                                             std::vector<std::vector<std::vector<double>>>& F)
{
    std::vector<double> oldF;
    F.resize(2);
    F[0].resize(8);
    F[1].resize(8);
    for (int i = 0; i < 8; i++)
    {
        F[0][i].resize(dyn[0].z.size());
        F[1][i].resize(dyn[0].z.size());
    };
    if (step == 0)
    {
        for (int i = 0; i < 2; i++)
            LinacTools::calculateF(dyn[i], alpha, F[i]);
    }
    else
        LinacTools::calculateF(dyn[1], alpha, F[1]);

    // dyn[0] = dyn[1];

    for (int i = 0; i < dyn[0].z.size(); i++)
    {

        if (dyn[1].cellNumber[i] == -1)
            continue;

        dyn[1].p[i] = dyn[1].p[i] + (dt / 2.0) * (3 * F[1][0][i] - F[0][0][i]);
        dyn[1].z[i] = dyn[1].z[i] + (dt / 2.0) * (3 * F[1][1][i] - F[0][1][i]);

        if (dyn[1].z[i] < 0)
            continue;

        dyn[1].S11x[i] = dyn[1].S11x[i] + (dt / 2.0) * (3 * F[1][2][i] - F[0][2][i]);
        dyn[1].S12x[i] = dyn[1].S12x[i] + (dt / 2.0) * (3 * F[1][3][i] - F[0][3][i]);
        dyn[1].S22x[i] = dyn[1].S22x[i] + (dt / 2.0) * (3 * F[1][4][i] - F[0][4][i]);

        dyn[1].S11y[i] = dyn[1].S11y[i] + (dt / 2.0) * (3 * F[1][5][i] - F[0][5][i]);
        dyn[1].S12y[i] = dyn[1].S12y[i] + (dt / 2.0) * (3 * F[1][6][i] - F[0][6][i]);
        dyn[1].S22y[i] = dyn[1].S22y[i] + (dt / 2.0) * (3 * F[1][7][i] - F[0][7][i]);

        /*if (dyn[1].S11x[i] < 1e-9 && dyn[1].S22x[i]>0)
                dyn[1].S22x[i] = -dyn[1].S22x[i];

        if (dyn[1].S11y[i] < 1e-9 && dyn[1].S22y[i]>0)
                dyn[1].S22y[i] = -dyn[1].S22y[i];*/

        /*	if (dyn[1].S11x[i] < 0 && dyn[1].S22x[i]>0)
                {
                        dyn[1].S11x[i] = -dyn[1].S11x[i];
                        dyn[1].S22x[i] = -dyn[1].S22x[i];
                }


                if (dyn[1].S11y[i] < 0 && dyn[1].S22y[i]>0)
                {
                        dyn[1].S11y[i] = -dyn[1].S11y[i];
                        dyn[1].S22y[i] = -dyn[1].S22y[i];
                }*/
    };
    F[0] = F[1];
};

void LinacTools::updateMomentumsAndPositions(std::vector<std::vector<double>>& dyn,
                                             std::vector<double>& oneParticleData, double dt,
                                             double alpha, int step)
{
    // dyn[0] = dyn[1];

    if (step == 0)
    {
        for (int i = 0; i < 2; i++)
            LinacTools::calculateF(dyn[i], oneParticleData, alpha, FnOnePartilce[i]);
    }
    else
        LinacTools::calculateF(dyn[1], oneParticleData, alpha, FnOnePartilce[1]);

    for (int i    = 0; i < 8; i++)
        dyn[1][i] = dyn[1][i] + (dt / 2.0) * (3 * FnOnePartilce[1][i] - FnOnePartilce[0][i]);

    FnOnePartilce[0] = FnOnePartilce[1];
};

void LinacTools::updateMomentumsAndPositionsCorrect(std::vector<std::vector<double>>& dyn,
                                                    std::vector<double>& oneParticleData, double dt,
                                                    double alpha, int step)
{
    LinacTools::calculateF(dyn[1], oneParticleData, alpha, FnOnePartilce[1]);

    for (int i    = 0; i < 8; i++)
        dyn[1][i] = dyn[0][i] + (dt / 2.0) * (FnOnePartilce[1][i] + FnOnePartilce[0][i]);

    FnOnePartilce[0] = FnOnePartilce[1];
};

double LinacTools::InitialFitnessFunction(const std::vector<double>&      population,
                                          LinacTools::RFQFitnessDynamics* dynObj)
{
    /*std::vector<std::vector<std::vector<double>>> RFQControls;
    LinacTools::GenerateRFQControlsbyParameterss(RFQControls, population);

    int nCellRegular = population[0] + population[3] + population[6] - 2;
    int nCellMatcher = dynObj->RFQParameters[8];
    int nCell = nCellRegular + nCellMatcher;


    int succes;

    dynObj->outputData.clear();
    LinacTools::GenerateRFQForFlow(succes, dynObj->RFQParameters, dynObj->flows[0], RFQControls,
    dynObj->CellsLengths,
    dynObj->MinumalRadii, dynObj->u1c, dynObj->u2c, dynObj->MatcherRadii, dynObj->AllRadii,
    dynObj->AvEnergies);
    dynObj->dynamicsSimulator.SimulateRFQFlowDynamics(0, nCell, dynObj->flows[0],
    dynObj->RFQParameters, dynObj->u1c,
    dynObj->u2c, dynObj->CellsLengths, dynObj->AllRadii, dynObj->outputData, RFQControls[1][1][0]);
    return
    dynObj->outputData[0]->dataAdd[3].back();*/
    return 0;
};

void LinacTools::GenerateRFQControlsbyParameterss(
    const std::vector<double>& RFQParameters, std::vector<double>& Modulations,
    std::vector<double>& SyncPhases, std::vector<double>& MatcherRadii,
    std::vector<double>& OutputRadii, std::vector<std::vector<double>> RFQApproxParameters)
{

    //"Number of cells", "Output phase", "Output Modulation", "Approximation degree"

    double Phase0      = -commtools::PI() / 2;
    double Modulation0 = 1;
    double Phase1;
    double Modulation1;
    int    ncells;
    double degree;
    double degree1;
    Modulations.clear();
    SyncPhases.clear();
    OutputRadii.clear();
    int nCellMatcher    = int(RFQApproxParameters[4][0]);
    int nCellMatcherOut = int(RFQApproxParameters[5][0]);

    for (int i = 0; i < nCellMatcher; i++)
    {
        SyncPhases.push_back(-commtools::PI() / 2);
        Modulations.push_back(1);
    };
    for (int section = 0; section < RFQApproxParameters.size() - 2; section++)
    {
        Phase1        = RFQApproxParameters[section][1];
        Modulation1   = RFQApproxParameters[section][2];
        ncells        = int(RFQApproxParameters[section][0]);
        degree        = RFQApproxParameters[section][3];
        degree1       = RFQApproxParameters[section][4];
        double kMod   = (Modulation1 - Modulation0) / pow(double(ncells), degree1);
        double kPhase = (Phase1 - Phase0) / pow(double(ncells), degree);
        for (int i = 0; i < ncells; i++)
        {

            Modulations.push_back(Modulation0 + pow(double(i + 1), degree1) * kMod);
            SyncPhases.push_back(Phase0 + pow(double(i + 1), degree) * kPhase);
        };
        Modulation0 = Modulation1;
        Phase0      = Phase1;
    }

    double mD = 2;

    double Rmax = RFQApproxParameters[4][1];

    double kMatcher = (Rmax - RFQParameters[1]) / (pow(nCellMatcher, mD));

    MatcherRadii.resize(nCellMatcher);

    for (int i = 0; i < nCellMatcher; i++)
    {
        MatcherRadii[i] = RFQParameters[1] + pow(nCellMatcher - i, mD) * kMatcher;
    }

    double lph = SyncPhases.back();

    OutputRadii.resize(nCellMatcherOut);

    double RmaxOut = RFQApproxParameters[5][1];

    double kOut = (RmaxOut - RFQParameters[1]) / (pow(nCellMatcherOut, mD));

    for (int i = 0; i < nCellMatcherOut; i++)
    {
        SyncPhases.push_back(lph);
        Modulations.push_back(1);
        OutputRadii[i] = RFQParameters[1] + pow(i, mD) * kOut;
    };

    SyncPhases.push_back(SyncPhases.back());
};

void LinacTools::changeRFQcontrols(std::vector<std::vector<std::vector<double>>>& RFQControlsInput)
{
    /*int ncell = RFQControlsInput[0][0].back();

    ncell = (ncell - 50) + rand() % (ncell + 100);

    int dn1 = ncell / (RFQControlsInput[0][0].size() - 1);*/

    for (int i = 0; i < RFQControlsInput.size(); i++)
    {

        /*for (int ii = 0; ii < RFQControlsInput[i][0].size() - 2; ii++)
                RFQControlsInput[i][0][ii + 1] = RFQControlsInput[i][0][ii] + dn1;

        RFQControlsInput[i][0].back() = ncell;*/

        for (int j = 1; j < RFQControlsInput[i][1].size() - 2; j++)
        {
            int r = rand() % 100;
            // distribution_lt.mean = RFQControlsInput[i][1][j];
            double mean = RFQControlsInput[i][1][j];
            if (r > 50)
            {
                double sigma = std::abs((RFQControlsInput[i][1][j + 1] - RFQControlsInput[i][1][j]) / 3);
                double dV    = std::abs(sigma * distribution_lt(generator_lt) + mean);
                RFQControlsInput[i][1][j] + dV;
            }
            else
            {
                double sigma = std::abs((RFQControlsInput[i][1][j] - RFQControlsInput[i][1][j - 1]) / 3);
                double dV    = std::abs(sigma * distribution_lt(generator_lt) + mean);
                RFQControlsInput[i][1][j] - dV;
            }
        };

        /*for (int j = 1; j < RFQControlsInput[i][1].size() - 1; j++)
        {

                double sigmaX = (RFQControlsInput[i][0][j + 1] - RFQControlsInput[i][0][j]) / 3;
                double meanX = RFQControlsInput[i][0][j];
                RFQControlsInput[i][0][j] = sigmaX * distribution_lt(generator_lt) + meanX;
        }*/
    }
};

void LinacTools::updateMomentums(LinacDynamicsLong& dyn, double dt, double alpha)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        if (dyn.cellNumber[i] == -1)
            continue;

        dyn.p[i] = dyn.p[i] + dt * alpha * dyn.Ez[i];
    };
};
void LinacTools::updatePositions(LinacDynamicsLong& dyn, double dt, double alpha)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        // if (dyn.cellNumber[i] == -1)
        //	continue;

        dyn.z[i] = dyn.z[i] + dt * dyn.p[i];

        if (dyn.z[i] < 0)
            continue;
    }
};

void LinacTools::InCell(const std::vector<double>& z, LinacDynamicsTransv& dyn,
                        std::vector<double>& W1, std::vector<double>& W2, std::vector<double>& rho)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        int cell = dyn.cellNumber[i];
        if (cell == -1)
            continue;
        if (dyn.cellNumber[i] == -2)
        {
            for (int j = 0; j < z.size() - 1; j++)
            {
                if (dyn.z[i] >= z[j] && dyn.z[i] <= z[j + 1])
                {
                    dyn.cellNumber[i] = j;
                    double dz         = z[j + 1] - z[j];
                    W1[i]             = (z[j + 1] - dyn.z[i]) / dz;
                    break;
                }
            };
            continue;
        }

        int lastCell = dyn.cellNumber[i];

        if (dyn.z[i] >= z[lastCell] && dyn.z[i] <= z[lastCell + 1])
        {
            dyn.cellNumber[i] = lastCell;

            double dz = z[lastCell + 1] - z[lastCell];
            W1[i]     = (z[lastCell + 1] - dyn.z[i]) / dz;

            continue;
        };

        dyn.cellNumber[i] = -1;

        for (int k = 1; k < 3; k++)
        {
            if (dyn.z[i] >= z[lastCell + k] && dyn.z[i] <= z[lastCell + k + 1])
            {
                dyn.cellNumber[i] = lastCell + k;

                double dz = z[lastCell + 1 + k] - z[lastCell + k];
                W1[i]     = (z[lastCell + 1 + k] - dyn.z[i]) / dz;

                break;
            };
        }

        if (dyn.cellNumber[i] == -1)
        {
            for (int k = -1; k > -4; k--)
            {
                if (lastCell + k < 0)
                    break;
                if (dyn.z[i] >= z[lastCell + k] && dyn.z[i] <= z[lastCell + k + 1])
                {
                    dyn.cellNumber[i] = lastCell + k;

                    double dz = z[lastCell + 1 + k] - z[lastCell + k];
                    W1[i]     = (z[lastCell + 1 + k] - dyn.z[i]) / dz;
                    break;
                };
            }
        };
    };
};

void LinacTools::InCell(const std::vector<double>& z, std::vector<double>& dyn,
                        std::vector<double>& oneParticleData, RFQExternalGrid& grid, double time,
                        double freq, double ph0, double I)
{
    int    cell = oneParticleData[0];
    double W1;
    double cost = std::cos(time * freq + ph0);

    if (cell == -1)
        return;

    double v  = dyn[0] * commtools::LIGHT_VELOCITY();
    double rx = sqrt(dyn[2]);
    double ry = sqrt(dyn[5]);

    double l = 5e-7;
    if (dyn[2] < l || dyn[5] < l)
    {
        rx = sqrt(l);
        ry = sqrt(l);
    }

    double KVForceX = I * (1 - (rx - ry) / (rx + ry)) /
                      (2 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * v * rx * ry);
    double KVForceY = I * (1 + (rx - ry) / (rx + ry)) /
                      (2 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * v * rx * ry);

    if (std::isinf(KVForceX) || std::isnan(KVForceX) || std::isinf(KVForceY) || std::isnan(KVForceY))
    {
        int tt = 0;
    }

    // KVForceX = 0;
    // KVForceY = 0;

    if (oneParticleData[0] == -2)
    {
        for (int j = 0; j < z.size() - 1; j++)
        {
            if (dyn[1] >= z[j] && dyn[1] <= z[j + 1])
            {
                oneParticleData[0] = j;
                double dz          = z[j + 1] - z[j];
                W1                 = (z[j + 1] - dyn[1]) / dz;

                oneParticleData[1] =
                    grid.Ez[oneParticleData[0]] * W1 + grid.Ez[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[1] = oneParticleData[1] * cost;

                oneParticleData[2] =
                    grid.Ex[oneParticleData[0]] * W1 + grid.Ex[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[2] = oneParticleData[2] * cost + KVForceX;

                oneParticleData[3] =
                    grid.Ey[oneParticleData[0]] * W1 + grid.Ey[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[3] = oneParticleData[3] * cost + KVForceY;

                return;
            }
        };
        return;
    }

    int lastCell = oneParticleData[0];

    if (dyn[1] >= z[lastCell] && dyn[1] <= z[lastCell + 1])
    {
        oneParticleData[0] = lastCell;

        double dz = z[lastCell + 1] - z[lastCell];
        W1        = (z[lastCell + 1] - dyn[1]) / dz;
        oneParticleData[1] =
            grid.Ez[oneParticleData[0]] * W1 + grid.Ez[oneParticleData[0] + 1] * (1 - W1);
        oneParticleData[1] = oneParticleData[1] * cost;

        oneParticleData[2] =
            grid.Ex[oneParticleData[0]] * W1 + grid.Ex[oneParticleData[0] + 1] * (1 - W1);
        oneParticleData[2] = oneParticleData[2] * cost + KVForceX;

        oneParticleData[3] =
            grid.Ey[oneParticleData[0]] * W1 + grid.Ey[oneParticleData[0] + 1] * (1 - W1);
        oneParticleData[3] = oneParticleData[3] * cost + KVForceY;

        return;
    };

    oneParticleData[0] = -1;

    for (int k = 1; k < 3; k++)
    {
        if (dyn[1] >= z[lastCell + k] && dyn[1] <= z[lastCell + k + 1])
        {
            oneParticleData[0] = lastCell + k;

            double dz = z[lastCell + 1 + k] - z[lastCell + k];
            W1        = (z[lastCell + 1 + k] - dyn[1]) / dz;
            oneParticleData[1] =
                grid.Ez[oneParticleData[0]] * W1 + grid.Ez[oneParticleData[0] + 1] * (1 - W1);
            oneParticleData[1] = oneParticleData[1] * cost;

            oneParticleData[2] =
                grid.Ex[oneParticleData[0]] * W1 + grid.Ex[oneParticleData[0] + 1] * (1 - W1);
            oneParticleData[2] = oneParticleData[2] * cost + KVForceX;

            oneParticleData[3] =
                grid.Ey[oneParticleData[0]] * W1 + grid.Ey[oneParticleData[0] + 1] * (1 - W1);
            oneParticleData[3] = oneParticleData[3] * cost + KVForceY;

            return;
        };
    }

    if (oneParticleData[0] == -1)
    {
        for (int k = -1; k > -4; k--)
        {
            if (lastCell + k < 0)
                break;
            if (dyn[1] >= z[lastCell + k] && dyn[1] <= z[lastCell + k + 1])
            {
                oneParticleData[0] = lastCell + k;

                double dz = z[lastCell + 1 + k] - z[lastCell + k];
                W1        = (z[lastCell + 1 + k] - dyn[1]) / dz;
                oneParticleData[1] =
                    grid.Ez[oneParticleData[0]] * W1 + grid.Ez[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[1] = oneParticleData[1] * cost;

                oneParticleData[2] =
                    grid.Ex[oneParticleData[0]] * W1 + grid.Ex[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[2] = oneParticleData[2] * cost + KVForceX;

                oneParticleData[3] =
                    grid.Ey[oneParticleData[0]] * W1 + grid.Ey[oneParticleData[0] + 1] * (1 - W1);
                oneParticleData[3] = oneParticleData[3] * cost + KVForceY;

                return;
            };
        }
    };
};

void LinacTools::FieldInterpolate(LinacDynamicsTransv& dyn, RFQExternalGrid& grid,
                                  std::vector<double>& W1, double time, double freq, double ph0,
                                  double I, int flag)
{
    double cost = std::cos(time * freq + ph0);

    for (int i = 0; i < dyn.z.size(); i++)
    {
        if (dyn.cellNumber[i] == -1)
            continue;

        double v  = dyn.p[i] * commtools::LIGHT_VELOCITY();
        double rx = sqrt(dyn.S11x[i]);
        double ry = sqrt(dyn.S11y[i]);

        double l = 1e-8;

        if (dyn.S11x[i] < l)
            rx = sqrt(l);

        if (dyn.S11y[i] < l)
            ry = sqrt(l);

        double KVForceX = I * (1 - (rx - ry) / (rx + ry)) /
                          (2 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * v * rx * ry);
        double KVForceY = I * (1 + (rx - ry) / (rx + ry)) /
                          (2 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * v * rx * ry);

        dyn.Ex[i] =
            grid.Ex[dyn.cellNumber[i]] * W1[i] + grid.Ex[dyn.cellNumber[i] + 1] * (1 - W1[i]);
        dyn.Ex[i] = dyn.Ex[i] * cost + KVForceX;

        dyn.Ey[i] =
            grid.Ey[dyn.cellNumber[i]] * W1[i] + grid.Ey[dyn.cellNumber[i] + 1] * (1 - W1[i]);
        dyn.Ey[i] = dyn.Ey[i] * cost + KVForceY;

        dyn.Ez[i] =
            grid.Ez[dyn.cellNumber[i]] * W1[i] + grid.Ez[dyn.cellNumber[i] + 1] * (1 - W1[i]);

        if (flag == 1)
            dyn.EzCol[i] = grid.EzCol[dyn.cellNumber[i]] * W1[i] +
                           grid.EzCol[dyn.cellNumber[i] + 1] * (1 - W1[i]);

        dyn.Ez[i] = dyn.Ez[i] * cost + dyn.EzCol[i];
    }
};

void LinacTools::InCell(const std::vector<double>& z, LinacDynamicsLong& dyn,
                        std::vector<double>& W1, std::vector<double>& W2, std::vector<double>& rho)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        int cell = dyn.cellNumber[i];
        if (cell == -1)
            continue;
        if (dyn.cellNumber[i] == -2)
        {
            for (int j = 0; j < z.size() - 1; j++)
            {
                if (dyn.z[i] >= z[j] && dyn.z[i] <= z[j + 1])
                {
                    dyn.cellNumber[i] = j;
                    double dz         = z[j + 1] - z[j];
                    W1[i]             = (z[j + 1] - dyn.z[i]) / dz;

                    rho[j]     = rho[j] + dyn.longQ * W1[i];
                    rho[j + 1] = rho[j + 1] + dyn.longQ * (1 - W1[i]);
                    break;
                }
            };
            continue;
        }

        int lastCell = dyn.cellNumber[i];

        if (dyn.z[i] >= z[lastCell] && dyn.z[i] <= z[lastCell + 1])
        {
            dyn.cellNumber[i] = lastCell;

            double dz = z[lastCell + 1] - z[lastCell];
            W1[i]     = (z[lastCell + 1] - dyn.z[i]) / dz;

            rho[lastCell]     = rho[lastCell] + dyn.longQ * W1[i];
            rho[lastCell + 1] = rho[lastCell + 1] + dyn.longQ * (1 - W1[i]);

            continue;
        };

        dyn.cellNumber[i] = -1;

        for (int k = 1; k < 3; k++)
        {
            if (dyn.z[i] >= z[lastCell + k] && dyn.z[i] <= z[lastCell + k + 1])
            {
                dyn.cellNumber[i] = lastCell + k;

                double dz = z[lastCell + 1 + k] - z[lastCell + k];
                W1[i]     = (z[lastCell + 1 + k] - dyn.z[i]) / dz;

                rho[lastCell + k]     = rho[lastCell + k] + dyn.longQ * W1[i];
                rho[lastCell + 1 + k] = rho[lastCell + 1 + k] + dyn.longQ * (1 - W1[i]);

                break;
            };
        }

        if (dyn.cellNumber[i] == -1)
        {
            for (int k = -1; k > -4; k--)
            {
                if (lastCell + k < 0)
                    break;
                if (dyn.z[i] >= z[lastCell + k] && dyn.z[i] <= z[lastCell + k + 1])
                {
                    dyn.cellNumber[i] = lastCell + k;

                    double dz = z[lastCell + 1 + k] - z[lastCell + k];
                    W1[i]     = (z[lastCell + 1 + k] - dyn.z[i]) / dz;

                    rho[lastCell + k]     = rho[lastCell + k] + dyn.longQ * W1[i];
                    rho[lastCell + 1 + k] = rho[lastCell + 1 + k] + dyn.longQ * (1 - W1[i]);

                    break;
                };
            }
        };
    };
};

void LinacTools::FieldInterpolate(LinacDynamicsLong& dyn, RFQExternalGrid& grid,
                                  std::vector<double>& W1, double time, double freq, double ph0,
                                  double I)
{
    double cost = std::cos(time * freq + ph0);
    for (int i = 0; i < dyn.z.size(); i++)
    {
        if (dyn.cellNumber[i] == -1)
            continue;

        dyn.Ez[i] =
            grid.Ez[dyn.cellNumber[i]] * W1[i] + grid.Ez[dyn.cellNumber[i] + 1] * (1 - W1[i]);
        dyn.Ez[i] = dyn.Ez[i] * cost;
    }
};

double f2(double m, double k, double r, double t)
{
    return t -
           (commtools::PI() / 4) * (m * m - 1) /
               (m * m * boost::math::cyl_bessel_i(0, 2 * k * r / (m + 1)) +
                boost::math::cyl_bessel_i(0, 2 * m * k * r / (m + 1)));
}

void LinacTools::LinacDynSimulator::SetFieldSolverParameters()
{
    std::vector<double> coloumbSolverParameters(12);

    coloumbSolverParameters[0] = 1.9;
    coloumbSolverParameters[1] = 0.001;
    coloumbSolverParameters[2] = 1.7;
    coloumbSolverParameters[3] = 0.01;

    coloumbSolverParameters[4] = 1e-7;
    coloumbSolverParameters[5] = 1e-6;

    coloumbSolverParameters[6] = 5;
    coloumbSolverParameters[7] = parameters[5];

    coloumbSolverParameters[8]  = 0;
    coloumbSolverParameters[9]  = 1;
    coloumbSolverParameters[10] = 1;
    coloumbSolverParameters[11] = 1;

    //	coloumbSolver.SetParameters(coloumbSolverParameters);
};

LinacTools::LinacDynSimulator::LinacDynSimulator()
{
    parameters.resize(7);
    parameters[0] = 50;
    parameters[1] = 50;
    parameters[2] = 20;
    parameters[3] = 2;
    parameters[4] = 0.004;
    parameters[5] = 2;
    parameters[6] = 0;
    SetFieldSolverParameters();
};

void LinacTools::LinacDynSimulator::Init(){

};

/*
void LinacTools::RFQExternalGrid::CreateRFQGrid(const std::vector<double>& u1cells,
                                                const std::vector<double>& L,
                                                const std::vector<double>& LL,
                                                const std::vector<double>& MinumalRadii, int nCells,
                                                int Ltotal, double lambda, double dz,
                                                double voltage)
{

    std::vector<double> LL(nCells + 1);

    LL[0] = 0;

    for (int i    = 0; i < nCells; i++)
        LL[i + 1] = LL[i] + device->GetSequencesOfParametersCalc()[0][i];

    double Ltotal = LL.back();

    z.clear();
    Ez.clear();
    Ex.clear();
    Ey.clear();
    double dzGrid;
    double zCurrent = -6 * lambda * L[0];
    int    period   = 0;
    double eUTmax   = voltage;

    while (1)
    {
        z.push_back(zCurrent);

        double kk;
        if (period % 2 == 0)
            kk = commtools::PI();
        else
            kk = 0;

        double EzCurrent;

        double K = commtools::PI() / L[period];

        double A = (2 * u1cells[period] * eUTmax / L[period]);

        EzCurrent = A * std::cos(K * (zCurrent - LL[period]) + kk);

        double k = 1 - 4 * u1cells[period] / commtools::PI();

        double B = eUTmax * k / (MinumalRadii[period] * MinumalRadii[period]);

        double C = commtools::PI() * eUTmax * u1cells[period] / (L[period] * L[period]);

        // C = 0;

        double ExCurrent = B + C * std::sin(K * (zCurrent - LL[period]) + kk);
        double EyCurrent = -B + C * std::sin(K * (zCurrent - LL[period]) + kk);

        // ExCurrent = 0;
        // EyCurrent = 0;
        if (zCurrent < 0)
        {
            EzCurrent = 0;
            ExCurrent = 0;
            EyCurrent = 0;
        }

        Ez.push_back(EzCurrent);
        Ex.push_back(ExCurrent);
        Ey.push_back(EyCurrent);

        dzGrid = L[period] / dz;

        zCurrent = zCurrent + dzGrid;

        if (zCurrent > LL[period + 1])
            period++;

        if (period >= nCells)
            break;
    };

    // empty space after electrodes
    while (1)
    {
        z.push_back(zCurrent);
        Ez.push_back(0);
        Ex.push_back(0);
        Ey.push_back(0);

        zCurrent = zCurrent + dzGrid;
        if (zCurrent > LL.back() + 2 * L.back())
            break;
    };
    EzCol.resize(z.size());
    rho.resize(z.size());
};*/

void LinacTools::LinacDynSimulator::CalculateRFQAcceptance(
    int flagType, LinacFlow& flow, const std::vector<double>& RFQParameters,
    const std::vector<double>& L, const std::vector<double>& Modulations,
    const std::vector<double>& MinumalRadii, const std::vector<double>& MatcherRadii,
    const std::vector<double>& AccEff, std::vector<DynamicsData*>& outputData)
{
    int np1               = flow.nParticlesLong;
    int np2               = flow.nParticlesTransv;
    flow.nParticlesLong   = 1;
    flow.nParticlesTransv = 1;

    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCells       = nCellRegular + nCellMatcher;

    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];
    double dt     = lambda / parameters[0];
    double dz     = 0.5 * parameters[0];

    std::vector<double> LL(nCells + 1);

    LL[0] = 0;

    for (int i    = 0; i < nCells; i++)
        LL[i + 1] = LL[i] + L[i];

    double Ltotal = LL.back();

    externalGrid.CreateRFQGrid(AccEff, L, LL, MinumalRadii, nCells, lambda, dz, flow.voltage);
    W1.resize(flow.nParticlesLong);
    W2.resize(flow.nParticlesLong);

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;

    std::vector<double> em;
    std::vector<double> ph0V;

    int npH = 1;

    double dph = commtools::PI() / 4;

    FnOnePartilce.resize(2);
    FnOnePartilce[0].resize(8);
    FnOnePartilce[1].resize(8);

    int v2;

    for (int phI = 0; phI < npH; phI++)
    {
        double ph0 = 0 + phI * dph;

        ph0V.push_back(ph0);

        double X1  = 5e-4;
        double dX1 = 5e-4;

        double X2  = 5e-4;
        double dX2 = 5e-4;

        while (1)
        {
            outputData.clear();

            flow.Xmax       = X2;
            flow.dXmax      = dX1;
            flow.emittanceX = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();
            flow.emittanceY = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();

            SimulateOneParticleRFQ(1, nCellMatcher, 1, nCells, flow, RFQParameters, AccEff, L, LL,
                                   MinumalRadii, v2, ph0, outputData);
            //	SimulateRFQFlowDynamics(nCellMatcher, 1, nCells, flow, RFQParameters, u1cells,
            // u2cells, L, MinumalRadii,
            // outputData, ph0);

            //	return;

            //	v2 = int(outputData[0]->dataAdd[2].back());

            if (std::abs(v2) < 0.001)
                break;
            else
                X2 = X2 * 2;
        };

        while (1)
        {
            outputData.clear();

            flow.Xmax       = (X1 + X2) / 2;
            flow.dXmax      = dX1;
            flow.emittanceX = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();
            flow.emittanceY = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();

            SimulateOneParticleRFQ(1, nCellMatcher, 1, nCells, flow, RFQParameters, AccEff, L, LL,
                                   MinumalRadii, v2, ph0, outputData);

            //	SimulateRFQFlowDynamics(nCellMatcher, 1, nCells, flow, RFQParameters, u1cells,
            // u2cells, L, MinumalRadii,
            // outputData, ph0);

            //	v2 = int(outputData[0]->dataAdd[2].back());

            if (v2 > 0)
                X1 = flow.Xmax;
            else
                X2 = flow.Xmax;

            if (std::abs(X1 - X2) < 1e-4)
            {
                flow.Xmax = X1;
                break;
            }
        }

        while (1)
        {

            outputData.clear();

            flow.dXmax      = dX2;
            flow.emittanceX = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();
            flow.emittanceY = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();

            SimulateOneParticleRFQ(1, nCellMatcher, 1, nCells, flow, RFQParameters, AccEff, L, LL,
                                   MinumalRadii, v2, ph0, outputData);

            //	SimulateRFQFlowDynamics(nCellMatcher, 1, nCells, flow, RFQParameters, u1cells,
            // u2cells, L, MinumalRadii,
            // outputData, ph0);

            //	v2 = int(outputData[0]->dataAdd[2].back());

            if (std::abs(v2) < 0.001)
                break;
            else
                dX2 = dX2 * 2;
        };

        while (1)
        {

            outputData.clear();

            flow.dXmax      = (dX1 + dX2) / 2;
            flow.emittanceX = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();
            flow.emittanceY = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();

            SimulateOneParticleRFQ(1, nCellMatcher, 1, nCells, flow, RFQParameters, AccEff, L, LL,
                                   MinumalRadii, v2, ph0, outputData);

            //	SimulateRFQFlowDynamics(nCellMatcher, 1, nCells, flow, RFQParameters, u1cells,
            // u2cells, L, MinumalRadii,
            // outputData, ph0);

            //	v2 = int(outputData[0]->dataAdd[2].back());

            if (v2 > 0)
                dX1 = flow.dXmax;
            else
                dX2 = flow.dXmax;

            if (std::abs(dX1 - dX2) < 1e-4)
            {
                flow.dXmax = dX1;
                break;
            }
        }

        flow.emittanceX = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();
        flow.emittanceY = 1e6 * flow.Xmax * flow.dXmax * flow.GetStartVelocity();

        flow.allParameters[3] = flow.emittanceX;
        flow.allParameters[4] = flow.emittanceX;
        em.push_back(flow.emittanceX);
    }
    //	flow.allParameters[];

    double rr             = flow.channelRelRadius;
    flow.channelRelRadius = 1;

    SimulateOneParticleRFQ(-1, nCellMatcher, 1, nCells, flow, RFQParameters, AccEff, L, LL,
                           MinumalRadii, v2, 0, outputData);
    flow.channelRelRadius = rr;
    double betaX          = flow.oneParticle[2] / flow.emittanceX;
    double AlphaX         = -flow.oneParticle[3] / flow.emittanceX;
    double gammaX         = flow.oneParticle[4] / flow.emittanceX;

    double betaY  = flow.oneParticle[5] / flow.emittanceY;
    double AlphaY = -flow.oneParticle[6] / flow.emittanceY;
    double gammaY = flow.oneParticle[7] / flow.emittanceY;

    // double Xmax = sqrt(betaX*flow.emittanceX);
    double Xmax = sqrt(flow.oneParticle[2]);

    // double dXmaxL = sqrt(gammaX*flow.emittanceX);

    double dXmaxL = sqrt(flow.oneParticle[4]);

    double dXmax = dXmaxL / flow.GetStartVelocity();

    flow.Xmax  = Xmax;
    flow.dXmax = dXmax;

    flow.Ymax  = Xmax;
    flow.dYmax = dXmax;

    flow.allParameters[10] = Xmax;
    flow.allParameters[11] = dXmax;

    flow.nParticlesLong   = np1;
    flow.nParticlesTransv = np2;
}

void LinacTools::EzFieldSimulate(int centralCell, int ncells, RFQExternalGrid& grid, double Rb,
                                 double Rch)
{
    if (centralCell == -1)
        return;
    dZtmp.resize(2 * ncells + 2);

    int nDz = 2 * ncells + 2;

    std::vector<double> V;

    double Vstart = 300;
    int    k      = 0;

    double dr = 1e-4;
    double dz = grid.z[centralCell] - grid.z[centralCell - 1];

    int nDr = Rch / dr;
    V.resize(nDz * nDr);
    Vtmp.resize(nDz * nDr);
    rightV.resize(nDz * nDr);

    for (int i = centralCell - ncells; i < centralCell + ncells + 2; i++)
    {
        double dZtmp = grid.z[i + 1] - grid.z[i];
        grid.rho[i]  = grid.rho[i] / dZtmp;
    }
    for (int j = 0; j < nDr; + j++)
    {
        for (int i = centralCell - ncells; i < centralCell + ncells + 2; i++)
        {
            rightV[k] = grid.rho[i] / commtools::VACUUM_PERMITTIVITY();

            if (j * dr > Rb)
                rightV[k] = 0;

            if (std::isnan(rightV[k]) || std::isinf(rightV[k]))
            {
                int tt = 0;
            };

            // rightV[k] = 0;
            k++;
        }
    }

    double tol   = 0;
    double omega = 1.7;
    int    iter  = 0;

    for (int i = 0; i < nDz; i++)
    {
        Vtmp[nDz * (nDr - 1) + i] = 0;
    }

    while (1)
    {
        tol = 0;
        for (int i = 0; i < nDz; i++)
        {
            Vtmp[i] = Vtmp[nDz + i];
            // Vtmp[i] = 0;
        }

        for (int j = 1; j < nDr - 1; j++)
        {
            Vtmp[nDz * j] = Vtmp[nDz * j + 1];
            // Vtmp[nDz*j] = 0;
            for (int i = 1; i < 2 * ncells + 1; i++)
            {

                int    ind = i + j * nDz;
                int    l   = i + (j - 1) * nDz;
                int    u   = i + (j + 1) * nDz;
                double r   = j * dr;
                double tmp = Vtmp[ind];

                double vv = rightV[ind] + (Vtmp[ind - 1] + Vtmp[ind + 1]) / (dz * dz) +
                            (Vtmp[l] + Vtmp[u]) / (dr * dr) + (Vtmp[u] - Vtmp[l]) / (2 * dr * r);

                double vvv = 2.0 / (dz * dz) + 2.0 / (dr * dr);

                Vtmp[ind] = Vtmp[ind] * (1 - omega) + omega * vv / vvv;

                if (tmp != 0 && std::abs(Vtmp[ind] - tmp) / tmp > tol)
                    tol = std::abs(Vtmp[ind] - tmp) / tmp;
            };

            Vtmp[nDz * (j + 1) - 1] = Vtmp[nDz * (j + 1) - 2];
            //	Vtmp[nDz*(j + 1) - 1] = 0;
        }

        // Vtmp[2 * ncells + 1] = VV[2 * ncells];
        iter++;
        if (iter > 5 && tol < 1e-3)
            break;
    };
    k = 1;
    for (int i = centralCell - ncells + 1; i < centralCell + ncells - 1; i++)
    {
        grid.EzCol[i] = -(Vtmp[k + 1] - Vtmp[k - 1]) / (2 * dz);
        if (std::isnan(grid.EzCol[i]) || std::isinf(grid.EzCol[i]) || grid.EzCol[i])
        {
            int tt = 0;
        };
        k++;
    };
};

void LinacTools::LinacDynSimulator::SimulateRFQFlowDynamics(
    int flagType, int nCells, LinacFlow& flow, const std::vector<double>& RFQParameters,
    const std::vector<double>& L, const std::vector<double>& MinumalRadii,
    const std::vector<double>& AccEff, std::vector<DynamicsData*>& outputData, double startPhase)
{

    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];
    double dt     = lambda / parameters[0];
    double dz     = 0.5 * parameters[0];

    std::vector<double> LL(nCells + 1);

    LL[0] = 0;

    for (int i    = 0; i < nCells; i++)
        LL[i + 1] = LL[i] + L[i];

    double Ltotal = LL.back();

    externalGrid.CreateRFQGrid(AccEff, L, LL, MinumalRadii, nCells, lambda, dz, flow.voltage);
    W1.resize(flow.nParticlesTransv);
    W2.resize(flow.nParticlesTransv);

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;
    double ph0   = startPhase;

    SimulateRFQFlowDynamics(0, flagType, nCells, flow, RFQParameters, L, MinumalRadii, AccEff,
                            outputData, startPhase);
};

void LinacTools::LinacDynSimulator::SimulateRFQFlowDynamics(
    int startCell, int flagType, int nCells, LinacFlow& flow,
    const std::vector<double>& RFQParameters, const std::vector<double>& L,
    const std::vector<double>& MinumalRadii, const std::vector<double>& AccEff,
    std::vector<DynamicsData*>& outputData, double startPhase)
{

    /*FILE* fp;
    fopen_s(&fp, "tt.txt", "w");
    double l = 1.031231465436546;
    fprintf(fp, "%.12lf", l);

    fclose(fp);*/

    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];
    double dt     = lambda / parameters[0];
    // startCell = 8;
    // startPhase = 0;
    double dz = 0.5 * parameters[0];

    double alpha = flow.GetAlpha();

    std::vector<double> LL(nCells + 1);

    LL[0] = 0;
    std::vector<std::vector<float>> OutParams(15);

    for (int i = 0; i < nCells; i++)
    {
        LL[i + 1] = LL[i] + L[i];
        OutParams[6].push_back(LL[i]);
        OutParams[7].push_back(MinumalRadii[i]);
    }

    flow.GenerateParticlesRFQ(lambda, LL[startCell], RFQParameters[1]);

    double Ltotal = LL.back();

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;

    int step = 0;

    DynamicsData* outputDataTmp;
    outputDataTmp = new DynamicsData();

    if (flagType == 1)
    {
        outputDataTmp->InitL(4, flow.mass, flow.charge, sizeof(double), 1, lambda, flow.emittanceX);

        std::vector<unsigned int> saveIndexes;
        for (int i = 0; i < flow.nParticlesTransv; i++)
            saveIndexes.push_back(i);
        outputDataTmp->AddSavedTraces(saveIndexes, 0);
    }

    std::vector<float> TimeOut;

    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    int    cell;
    double eUTmax = flow.voltage;

    std::vector<LinacDynamicsTransv> dyn(2);
    dyn[0] = flow.dynamicsTransv;
    dyn[1] = flow.dynamicsTransv;

    int flag = 0;
    int lC   = -10;

    while (1)
    {
        if (step % 1 == 0)
        {
            TimeOut.push_back(t / (1e-9 * commtools::LIGHT_VELOCITY()));

            // int n = flow.dynamicsLong.checkBadParticles(cell, OutParams, LL, MinumalRadii,
            // lambda);
            int n1 = dyn[1].checkBadParticles(cell, OutParams, LL, MinumalRadii, lambda,
                                              flow.channelRelRadius);

            if (lC != cell)
            {
                float d = 0;
                if (cell % 2)
                    d = commtools::PI();

                OutParams[14].push_back(t * wfreq + startPhase + d);
                lC = cell;
            }

            double beta = OutParams[4].back();

            // double defocusing = u1cells[cell] * eUTmax*std::abs(std::sin(u2cells[cell])) / (beta*beta);

            OutParams[5].push_back(0);

            //	if (n1 == 0 || dyn[1].cellNumber[0] ==-1)
            //		break;
            if (n1 == 0)
                break;
            if (flagType == 1)
            {
                outputDataTmp->SetData(dyn[1].GetData(), t / (commtools::LIGHT_VELOCITY()), 0, 1,
                                       1);
            }
        };

        // LinacTools::EzFieldSimulate(flow.dynamicsLong.cellNumber[0], dz, externalGrid);

        // int lC = dyn[1].cellNumber[0];
        LinacTools::InCell(externalGrid.z, dyn[1], W1, W2, externalGrid.rho);

        //	if (lC != dyn[1].cellNumber[0])
        //		OutParams[14].push_back(t / (commtools::LIGHT_VELOCITY()));

        if (step % 10 == 0 && flow.impulseCurrent != 0)
        {
            flag = 1;
            LinacTools::chargeCalculate(dyn[1], W1, externalGrid.rho, OutParams[1].back());
            LinacTools::EzFieldSimulate(dyn[1].cellNumber[0], dz, externalGrid, OutParams[1].back(),
                                        RFQParameters[1]);
            memset(&externalGrid.rho[0], 0, sizeof(externalGrid.rho[0]) * externalGrid.rho.size());
        };

        LinacTools::FieldInterpolate(dyn[1], externalGrid, W1, t, wfreq, startPhase,
                                     flow.impulseCurrent, flag);

        if (step % 10 == 0 && flow.impulseCurrent != 0)
            memset(&externalGrid.EzCol[0], 0,
                   sizeof(externalGrid.rho[0]) * externalGrid.rho.size());

        flag = 0;
        // LinacTools::updateMomentums(flow.dynamicsLong, dt, alpha);
        // LinacTools::updatePositions(flow.dynamicsLong, dt, alpha);

        if (step == 0)
            dyn[0] = dyn[1];

        LinacTools::updateMomentumsAndPositions(dyn, dt, alpha, step, F);

        t    = t + dt;
        step = step + 1;
    };

    outputDataTmp->SetAdditionalData(TimeOut, OutParams);
    outputData.push_back(outputDataTmp);

    /*if (flagType == 1)
    {
    }*/
}

void LinacTools::LinacDynSimulator::SimulateRFQLongDynamics(
    int startCell, int flagType, int nCells, LinacFlow& flow,
    const std::vector<double>& RFQParameters, const std::vector<double>& L,
    const std::vector<double>& MinumalRadii, const std::vector<double>& AccEff,
    std::vector<DynamicsData*>& outputData, double& Transmission, double startPhase)
{
    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];
    double dt     = lambda / parameters[0];
    double dz     = 0.5 * parameters[0];

    std::vector<double> LL(nCells + 1);

    LL[0] = 0;

    for (int i    = 0; i < nCells; i++)
        LL[i + 1] = LL[i] + L[i];

    double Ltotal = LL.back();

    externalGrid.CreateRFQGrid(AccEff, L, LL, MinumalRadii, nCells, lambda, dz, flow.voltage);
    W1.resize(flow.nParticlesLong);
    W2.resize(flow.nParticlesLong);

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;
    double ph0   = startPhase;

    flow.GenerateParticlesRFQ(lambda, LL[startCell], RFQParameters[1]);

    int                             step = 0;
    std::vector<std::vector<float>> OutParams(8);

    int    cell;
    double alpha = flow.GetAlpha();

    Transmission = 1;

    while (1)
    {
        if (step % 4 == 0)
        {

            int n = flow.dynamicsLong.checkBadParticles(cell, OutParams, LL, MinumalRadii, lambda);
            if (OutParams[3].back() < Transmission)
                Transmission = OutParams[3].back();

            if (n == 0)
                break;
        };

        LinacTools::InCell(externalGrid.z, flow.dynamicsLong, W1, W2, externalGrid.rho);
        LinacTools::FieldInterpolate(flow.dynamicsLong, externalGrid, W1, t, wfreq, startPhase,
                                     flow.impulseCurrent);

        LinacTools::updateMomentums(flow.dynamicsLong, dt, alpha);
        LinacTools::updatePositions(flow.dynamicsLong, dt, alpha);

        t    = t + dt;
        step = step + 1;
    };
}

void LinacTools::LinacDynSimulator::SimulateOneParticleRFQ(
    int timeSignum, int startCell, int flagType, int nCells, LinacFlow& flow,
    const std::vector<double>& RFQParameters, const std::vector<double>& u1cells,
    const std::vector<double>& L, const std::vector<double>& LL,
    const std::vector<double>& MinumalRadii, int& result, double startPhase,
    std::vector<DynamicsData*>& outputData)
{
    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];
    double dt     = timeSignum * lambda / parameters[0];
    double dz     = 0.5 * parameters[0];

    double alpha = flow.GetAlpha();

    flow.GenerateOneParticleRFQ(lambda, LL[startCell]);

    double Ltotal = LL.back();

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;

    int                              step = 0;
    std::vector<std::vector<double>> dyn(2);
    dyn[0] = flow.oneParticle;
    dyn[1] = flow.oneParticle;

    int lastCell = startCell;

    double rC;

    result = 1;

    /*DynamicsData* outputDataTmp;
    outputDataTmp = new DynamicsData();
    std::vector<float> TimeOut;
    std::vector<std::vector<float>> OutParams(6);*/

    while (1)
    {

        if (step % 4 == 0)
        {
            // TimeOut.push_back(t / (1e-9*commtools::LIGHT_VELOCITY()));

            for (int i = lastCell; i < LL.size(); i++)
            {

                if (LL[i] > dyn[1][1])
                {
                    lastCell = i;
                    break;
                };
            };

            if (lastCell >= LL.size() - 1)
                lastCell = LL.size() - 2;

            rC = std::max(dyn[1][2], dyn[1][5]);

            if (rC < 1e-8)
                rC = 1e-8;
            else
                rC = sqrt(rC);

            // OutParams[0].push_back(MinumalRadii[lastCell]);
            //	OutParams[1].push_back(rC);

            if (rC > MinumalRadii[lastCell] * flow.channelRelRadius || std::isinf(rC) || std::isnan(rC))
            {
                result           = 0;
                flow.oneParticle = dyn[1];
                break;
            };
        };

        LinacTools::InCell(externalGrid.z, dyn[1], flow.oneParticleData, externalGrid, t, wfreq,
                           startPhase, flow.impulseCurrent);

        LinacTools::updateMomentumsAndPositions(dyn, flow.oneParticleData, dt, alpha, step);

        // LinacTools::InCell(externalGrid.z, dyn[1], flow.oneParticleData, externalGrid, t, wfreq,
        // startPhase);

        //	LinacTools::updateMomentumsAndPositionsCorrect(dyn, flow.oneParticleData, dt, alpha,
        // step);

        t    = t + dt;
        step = step + 1;

        if (dyn[1][1] > Ltotal || dyn[1][1] < 0)
        {
            flow.oneParticle = dyn[1];
            break;
        }
    };
    flow.oneParticle = dyn[1];

    // outputDataTmp->SetAdditionalData(TimeOut, OutParams);
    // outputData.push_back(outputDataTmp);
}