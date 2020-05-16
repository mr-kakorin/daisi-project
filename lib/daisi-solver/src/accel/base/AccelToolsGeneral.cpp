#include "AccelToolsGeneral.h"
#include "AccelDeviceBase.h"
#include "../base/AccelFlow.h"
#include "Results.h"

#include <common_tools/constants.h>

std::vector<double>                           W1;
std::vector<double>                           W2;
std::vector<std::vector<std::vector<double>>> F;
void calculateF(LinacDynamicsTransv& dyn, double alpha, std::vector<std::vector<double>>& f)
{
    for (int i = 0; i < dyn.z.size(); i++)
    {
        double gamma = sqrt(1 + dyn.p[i] * dyn.p[i]);
        f[0][i]      = alpha * dyn.Ez[i];
        f[1][i]      = dyn.p[i] / gamma;
        f[2][i]      = 2 * dyn.S12x[i];
        f[3][i]      = dyn.Ex[i] * alpha * dyn.S11x[i] + dyn.S22x[i];
        f[4][i]      = 2 * dyn.Ex[i] * alpha * dyn.S12x[i];
        f[5][i]      = 2 * dyn.S12y[i];
        f[6][i]      = dyn.Ey[i] * alpha * dyn.S11y[i] + dyn.S22y[i];
        f[7][i]      = 2 * dyn.Ey[i] * alpha * dyn.S12y[i];
    }
}

void updateMomentumsAndPositions(std::vector<LinacDynamicsTransv>& dyn, double dt, double alpha,
                                 int step, std::vector<std::vector<std::vector<double>>>& F)
{
    std::vector<double> oldF;
    F.resize(2);
    F[0].resize(8);
    F[1].resize(8);
    for (int i = 0; i < 8; i++)
    {
        F[0][i].resize(dyn[0].z.size());
        F[1][i].resize(dyn[0].z.size());
    }

    if (step == 0)
    {
        for (int i = 0; i < 2; i++)
            calculateF(dyn[i], alpha, F[i]);
    }
    else
        calculateF(dyn[1], alpha, F[1]);

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
    };
    F[0] = F[1];
};

void FieldInterpolate(LinacDynamicsTransv& dyn, std::shared_ptr<ExternalAccelGrid>& grid,
                      std::vector<double>& W1, double time, double freq, double ph0, double I,
                      int flag)
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
            grid->Ex[dyn.cellNumber[i]] * W1[i] + grid->Ex[dyn.cellNumber[i] + 1] * (1 - W1[i]);
        dyn.Ex[i] = dyn.Ex[i] * cost + KVForceX;

        dyn.Ey[i] =
            grid->Ey[dyn.cellNumber[i]] * W1[i] + grid->Ey[dyn.cellNumber[i] + 1] * (1 - W1[i]);
        dyn.Ey[i] = dyn.Ey[i] * cost + KVForceY;

        dyn.Ez[i] =
            grid->Ez[dyn.cellNumber[i]] * W1[i] + grid->Ez[dyn.cellNumber[i] + 1] * (1 - W1[i]);

        if (flag == 1)
            dyn.EzCol[i] = grid->EzCol[dyn.cellNumber[i]] * W1[i] +
                           grid->EzCol[dyn.cellNumber[i] + 1] * (1 - W1[i]);

        dyn.Ez[i] = dyn.Ez[i] * cost + dyn.EzCol[i];
    }
};

void InCell(const std::vector<double>& z, LinacDynamicsTransv& dyn, std::vector<double>& W1,
            std::vector<double>& W2, std::vector<double>& rho)
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

void ExternalAccelGrid::CreateRFQGrid(const std::vector<double>& u1cells,
                                      const std::vector<double>& L,
                                      const std::vector<double>& MinumalRadii, int nCells,
                                      double lambda, double dz, double voltage)
{
    std::vector<double> LL(nCells + 1);

    LL[0] = 0;

    for (int i    = 0; i < nCells; i++)
        LL[i + 1] = LL[i] + L[i];

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
};
void ExternalAccelGrid::CreateDTLGrid(){

};

void SimulateFlowDynamicsEnvelopes(std::shared_ptr<ExternalAccelGrid> grid, int stepsPerCell,
                                   int startCell, int flagType, std::shared_ptr<AccelFlow>& flow,
                                   double frequency, const std::vector<double>& L,
                                   const std::vector<double>&           MinumalRadii,
                                   std::shared_ptr<SimulationDataAccel> outputData, int nsaveTr,
                                   double startPhase, double& progress)
{
    int nParticles = flow->GetnParticles();
    W1.resize(nParticles);
    W2.resize(nParticles);
    int    nParticlesSave = std::min(nParticles, nsaveTr);
    double lambda         = commtools::LIGHT_VELOCITY() / frequency;

    double dt = lambda / stepsPerCell;

    if (flagType == 1)
    {
        outputData->Init(8, flow->GetMass(), flow->GetCharge(), sizeof(double), nParticlesSave, 1,
                         0);
        outputData->AddBlock(0, nParticlesSave, 0, 1, 0);
        outputData->lambda     = lambda;
        outputData->emittanceX = flow->GetEmittanceX();
        outputData->emittanceY = flow->GetEmittanceY();
    };

    double alpha = flow->GetAlpha();

    std::vector<double> LL(MinumalRadii.size() + 1);

    LL[0] = 0;
    std::vector<std::vector<float>> OutParams(7);

    for (int i = 0; i < MinumalRadii.size(); i++)
    {
        LL[i + 1] = LL[i] + L[i];
        OutParams[5].push_back(LL[i]);
        OutParams[6].push_back(MinumalRadii[i]);
    }
    OutParams[5].push_back(LL.back());
    OutParams[6].push_back(MinumalRadii.back());

    flow->GenerateParticlesEnvelopes(lambda, LL[startCell], MinumalRadii.back());

    double Ltotal = LL.back();

    double t     = 0;
    double wfreq = 2 * commtools::PI() / lambda;

    int step = 0;

    std::vector<float> TimeOut;

    // channel R
    // beam R
    // Acceleration;
    // Transmission;
    int cell;

    std::vector<LinacDynamicsTransv> dyn(2);
    dyn[0] = flow->GetdynamicsTransv();
    dyn[1] = flow->GetdynamicsTransv();

    int flag = 0;
    int lC   = -10;

    while (1)
    {
        if (step % (stepsPerCell / 20) == 0)
        {
            TimeOut.push_back(t / (1e-9 * commtools::LIGHT_VELOCITY()));
            int n1 = dyn[1].checkBadParticles(cell, OutParams, LL, MinumalRadii, lambda,
                                              flow->GetchannelRelRadius());

            double beta = OutParams[4].back();

            if (n1 == 0)
                break;
            if (flagType == 1)
                outputData->SetData(dyn[1].GetData(), t / (commtools::LIGHT_VELOCITY()), 0, 1, 1);
        };
        progress = double(cell) / MinumalRadii.size();
        // LinacTools::EzFieldSimulaflow.dynamicsLong.cellNumber[0], dz, externalGrid);

        // int lC = dyn[1].cellNumber[0];
        InCell(grid->z, dyn[1], W1, W2, grid->rho);

        //	if (lC != dyn[1].cellNumber[0])
        //		OutParams[14].push_back(t / (commtools::LIGHT_VELOCITY()));

        /*	if (step % 10 == 0 && flow.impulseCurrent != 0)
                {
                        flag = 1;
                        LinacTools::chargeCalculate(dyn[1], W1, externalGrid.rho,
           OutParams[1].back());
                        LinacTools::EzFieldSimulate(dyn[1].cellNumber[0], dz, externalGrid,
           OutParams[1].back(),
           RFQParameters[1]); memset(&externalGrid.rho[0], 0,
           sizeof(externalGrid.rho[0])*externalGrid.rho.size());

                };*/

        FieldInterpolate(dyn[1], grid, W1, t, wfreq, startPhase, flow->GetimpulseCurrent(), flag);

        flag = 0;

        if (step == 0)
            dyn[0] = dyn[1];

        updateMomentumsAndPositions(dyn, dt, alpha, step, F);

        t    = t + dt;
        step = step + 1;
    };
    outputData->SetAdditionalData(TimeOut, OutParams);
    progress = 1;
};