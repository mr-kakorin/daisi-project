#include "DataTypes.h"
#include "Solver.h"
#include "armadillo"

#include <iostream>
#include <windows.h>
template void SolverPTI2daxsdouble::SimulateCPU<device2daxsdouble>(double& progress, bool& flagAbort,
                                                                   double             timeDilatation,
                                                                   device2daxsdouble& deviceStatus,
                                                                   std::vector<std::vector<DynamicsData*>>& outputData,
                                                                   std::mutex&                              plotMutex,
                                                                   SimulationData& simulationData, double outTime);

template void SolverPTI2daxsfloat::SimulateCPU<device2daxsfloat>(double& progress, bool& flagAbort,
                                                                 double timeDilatation, device2daxsfloat& deviceStatus,
                                                                 std::vector<std::vector<DynamicsData*>>& outputData,
                                                                 std::mutex& plotMutex, SimulationData& simulationData,
                                                                 double outTime);

template void SolverPTI2ddouble::SimulateCPU<device2ddouble>(double& progress, bool& flagAbort, double timeDilatation,
                                                             device2ddouble&                          deviceStatus,
                                                             std::vector<std::vector<DynamicsData*>>& outputData,
                                                             std::mutex& plotMutex, SimulationData& simulationData,
                                                             double outTime);

template void SolverPTI2dfloat::SimulateCPU<device2dfloat>(double& progress, bool& flagAbort, double timeDilatation,
                                                           device2dfloat&                           deviceStatus,
                                                           std::vector<std::vector<DynamicsData*>>& outputData,
                                                           std::mutex& plotMutex, SimulationData& simulationData,
                                                           double outTime);

template <class ParticleGridInterfaceType, class EmissionCurrentSolverType, class FieldSolver, class PointType>
template <class deviceType>
void SolverPTI<ParticleGridInterfaceType, EmissionCurrentSolverType, FieldSolver, PointType>::SimulateCPU(
    double& progress, bool& flagAbort, double timeDilatation, deviceType& deviceStatus,
    std::vector<std::vector<DynamicsData*>>& outputData, std::mutex& plotMutex, SimulationData& simulationData,
    double outTime)
{
    outTime = outTime * (1e-9) * commtools::LIGHT_VELOCITY();

    InitFieldSolver(deviceStatus);
    int memorySize = 1e5;

    simulationData.reset();

    std::vector<int> boundaryPoints = fieldSolver.getNearBoundarypoints_all();

    particleGridInterface.init(deviceStatus.GetGridData(), deviceStatus.mesh.templNumb, deviceStatus.mesh.flagMatrix,
                               boundaryPoints, deviceStatus.GetDomainBoundary(), memorySize);

    W.resize(deviceStatus.GetNumberParticlesFlows());
    for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
        W[i]   = new PointType[memorySize][9];

    Wtmp.resize(deviceStatus.GetNumberParticlesFlows());
    for (int i  = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
        Wtmp[i] = new PointType[memorySize][9];

    EmptyPlaces.resize(deviceStatus.GetNumberParticlesFlows());
    EmptyPlaces1.resize(deviceStatus.GetNumberParticlesFlows());

    std::vector<int> flag;

    emissionCurrentSolver.reset();

    for (int i = 0; i < deviceStatus.conductorList.size(); i++)
    {
        deviceStatus.conductorList[i].ResetPower();
    }
    int iteration = 0;

    std::vector<PointType> alpas(deviceStatus.GetNumberParticlesFlows());

    for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
    {
        flag.push_back(1);
        deviceStatus.GetFlow(i)->ReserveMemory(memorySize);
        emissionCurrentSolver.init(deviceStatus.GetFlow(i)->GetEmitterDevice()->GetParticleSource(),
                                   &particleGridInterface);
        deviceStatus.GetFlow(i)->GetDynamicsData()->clear();
        std::vector<unsigned int> v;
        deviceStatus.GetFlow(i)->GenerateParticles(v, 1, 1, 0);
        particleGridInterface.SearchStartCells(deviceStatus.GetFlow(i)->GetDynamicsData());
        alpas[i] = deviceStatus.GetFlow(i)->GetAlpha();
    }

    particleMover.init(alpas, 1);

    float maxEr;
    //	errors = estimateError(deviceStatus, maxEr);
    doIterartion(deviceStatus, plotMutex, simulationData, outputData, 1, iteration);

    if (emissionCurrentSolver.GetParameters() == 0)
        return;

    int              nFlows = deviceStatus.GetNumberParticlesFlows();
    std::vector<int> flowsSol;

    for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
    {
        if (deviceStatus.GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 1)
            nFlows--;
        else
            flowsSol.push_back(i);
    }

    std::vector<double> cathodeFields;
    std::vector<double> cathodeFieldsTmp1;

    int p1 = 0;
    int p2 = 7;

    for (int flowsIter = 0; flowsIter < nFlows; flowsIter++)
    {
        int i = flowsSol[flowsIter];
        if (deviceStatus.GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
        {
            cathodeFieldsTmp1 = emissionCurrentSolver.CalculateCathodeFields(
                deviceStatus.GetFlow(i)->GetEmitterDevice(), &particleGridInterface, deviceStatus.GetGridData(),
                particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY(), i, 0, deviceStatus.GetFlow(i)->GetMass(),
                deviceStatus.GetFlow(i)->GetCharge(), p2 - p1);
            for (int k = 0; k < cathodeFieldsTmp1.size(); k++)
            {
                cathodeFields.push_back(cathodeFieldsTmp1[k]);
            }
        }
    }

    std::vector<double> dp;
    double              dpp;
    deviceType          deviceStatusTmp;
    deviceType          deviceStatusTmp1;

    int dim = cathodeFields.size();

    std::vector<std::vector<double>> cathodeFieldsNew(dim);

    double*             A = new double[dim * dim];
    std::vector<double> b(dim);
    arma::mat           armaMat;
    arma::vec           armaVec;

    armaMat.ones(dim, dim);
    armaVec.ones(dim);

    double Emin = cathodeFields[0];

    for (int i = 0; i < dim; i++)
    {
        if (std::abs(cathodeFields[i]) < Emin)
            Emin = cathodeFields[i];
    }
    Emin = 0.0 * Emin;

    for (int i = 0; i < dim; i++)
    {
        armaVec(i) = (cathodeFields[i] - Emin) * (cathodeFields[i] - Emin);
    }

    deviceStatusTmp1 = deviceStatus;
    std::vector<double> ty;

    int kk = 0;

    std::vector<std::vector<double>> cathodeFieldsTmp;
    std::vector<std::vector<double>> FF;

    for (int flowsIter = 0; flowsIter < nFlows; flowsIter++)
    {
        int i = flowsSol[flowsIter];
        if (deviceStatus.GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
        {
            for (int j = p1; j < p2; j++)
            {
                //	deviceStatus.GetGridData()->densityReset();
                //	dpp = emissionCurrentSolver.ChangeCurrent(deviceStatusTmp.GetFlow(i)->GetEmitterDevice(), j,
                // 0.04*(1+pow(double(j),1.01)));
                double dpp1;

                dp.clear();
                int nDer = 2;
                cathodeFieldsTmp.resize(nDer);

                std::vector<std::vector<double>> F;
                F.resize(nDer + 1);

                for (int k = 0; k < cathodeFields.size(); k++)
                    F[0].push_back((cathodeFields[k] - Emin) * (cathodeFields[k] - Emin));

                float step      = 0.02;
                float stepStart = 2 * step;
                for (int Der = 0; Der < nDer; Der++)
                {
                    deviceStatus = deviceStatusTmp1;
                    dpp1         = emissionCurrentSolver.ChangeCurrent(deviceStatus.GetFlow(i)->GetEmitterDevice(), j,
                                                               step * (Der + 1));
                    dp.push_back(dpp1);
                    doIterartion(deviceStatus, plotMutex, simulationData, outputData, 0, iteration);

                    for (int i1flowsIter = 0; i1flowsIter < nFlows; i1flowsIter++)
                    {
                        int i1 = flowsSol[i1flowsIter];
                        // cathodeFieldsTmp =
                        // emissionCurrentSolver.CalculateCathodeFields(deviceStatusTmp.GetFlow(i1)->GetEmitterDevice(),
                        // &particleGridInterface, deviceStatusTmp.GetGridData(), particleMover.GetTimeStep(i1) /
                        // commtools::LIGHT_VELOCITY(), i1, 0, deviceStatusTmp.GetFlow(i1)->GetMass(),
                        // deviceStatusTmp.GetFlow(i1)->GetCharge(),p2-p1);
                        cathodeFieldsTmp1 = emissionCurrentSolver.CalculateCathodeFields(
                            deviceStatus.GetFlow(i1)->GetEmitterDevice(), &particleGridInterface,
                            deviceStatus.GetGridData(), particleMover.GetTimeStep(i1) / commtools::LIGHT_VELOCITY(), i1, 0,
                            deviceStatus.GetFlow(i1)->GetMass(), deviceStatus.GetFlow(i1)->GetCharge(), p2 - p1);

                        for (int k = 0; k < cathodeFieldsTmp1.size(); k++)
                        {
                            F[Der + 1].push_back((cathodeFieldsTmp1[k] - Emin) * (cathodeFieldsTmp1[k] - Emin));
                            //	cathodeFieldsTmp[Der].push_back(cathodeFieldsTmp1[k]);
                        }
                    }
                }

                while (1)
                {
                    int flagStop = 0;
                    for (int k = 0; k < p2 - p1; k++)
                    {
                        float par2 = std::abs(F[2][k + flowsIter * p2] - F[1][k + flowsIter * p2]);
                        float par1 = std::abs(F[1][k + flowsIter * p2] - F[0][k + flowsIter * p2]);
                        if (par2 > par1 && par2 - par1 > 0.05 * par2)
                            flagStop = 1;
                    };
                    if (flagStop == 0)
                        break;
                    F[1] = F[2];
                    step = step * 2;

                    deviceStatus = deviceStatusTmp1;
                    dpp1 =
                        emissionCurrentSolver.ChangeCurrent(deviceStatus.GetFlow(i)->GetEmitterDevice(), j, 2 * step);
                    dp.push_back(dpp1);
                    doIterartion(deviceStatus, plotMutex, simulationData, outputData, 0, iteration);

                    F[2].clear();
                    for (int i1flowsIter = 0; i1flowsIter < nFlows; i1flowsIter++)
                    {
                        int i1 = flowsSol[i1flowsIter];

                        // cathodeFieldsTmp =
                        // emissionCurrentSolver.CalculateCathodeFields(deviceStatusTmp.GetFlow(i1)->GetEmitterDevice(),
                        // &particleGridInterface, deviceStatusTmp.GetGridData(), particleMover.GetTimeStep(i1) /
                        // commtools::LIGHT_VELOCITY(), i1, 0, deviceStatusTmp.GetFlow(i1)->GetMass(),
                        // deviceStatusTmp.GetFlow(i1)->GetCharge(),p2-p1);
                        cathodeFieldsTmp1 = emissionCurrentSolver.CalculateCathodeFields(
                            deviceStatus.GetFlow(i1)->GetEmitterDevice(), &particleGridInterface,
                            deviceStatus.GetGridData(), particleMover.GetTimeStep(i1) / commtools::LIGHT_VELOCITY(), i1, 0,
                            deviceStatus.GetFlow(i1)->GetMass(), deviceStatus.GetFlow(i1)->GetCharge(), p2 - p1);

                        for (int k = 0; k < cathodeFieldsTmp1.size(); k++)
                            F[2].push_back((cathodeFieldsTmp1[k] - Emin) * (cathodeFieldsTmp1[k] - Emin));
                    }
                }

                for (int k = 0; k < F[0].size(); k++)
                {
                    float d = (-3 * F[0][k] + 4 * F[1][k] - F[2][k]) / (dpp1);
                    armaMat(k, kk) = d;
                }

                /*	FF = F;
                        for (int k = 0; k < FF[0].size(); k++)
                        {
                                FF[0][k] = 0.2*(3 * F[0][k] + 2 * F[1][k] + F[2][k] - F[4][k]);
                                FF[1][k] = 0.1*(4 * F[0][k] + 3 * F[1][k] + 2*F[2][k] + F[3][k]);
                                FF[2][k] = 0.2*(F[0][k] + F[1][k] + F[2][k] + F[3][k] + F[4][k]);
                                FF[3][k] = 0.1*(4*F[4][k] + 3*F[3][k] + 2*F[2][k] + F[1][k]);
                                FF[4][k] = 0.1*(3 * F[4][k] + 2 * F[3][k] + F[2][k] - F[0][k]);
                                armaMat(k, kk) = (-25.0*FF[0][k] + 48.0 * FF[1][k] - 36.0 * FF[2][k] + FF[3][k] - 3.0 *
                   FF[4][k]) / (12 * dp[0]); armaMat(k, kk) = (FF[4][k] - FF[0][k]) / (4 * dp[0]);
                        }*/
                //	A[dim*k + kk] = ((cathodeFieldsTmp1[k] - Emin) * (cathodeFieldsTmp1[k] - Emin) -
                //(cathodeFields[k] - Emin)* (cathodeFields[k] - Emin)) / (dpp1);

                //	armaMat(k, kk) = A[dim*k + kk];
                //	cathodeFieldsNew[kk].push_back(A[dim*k + kk]);

                kk++;
            }
        }
    }
    deviceStatus = deviceStatusTmp1;
    for (int i1 = 0; i1 < dim; i1++)
    {
        for (int i2 = 0; i2 < dim; i2++)
        {
            //			if (i1!=i2)
            //			armaMat(i1, i2) = 0;
        }
    }
    /*arma::mat armaMatE;
    armaMatE.zeros(dim, dim);


    for (int i1 = 0; i1 < dim; i1++)
            armaMatE(i1, i1) = 1;

    arma::mat armaMatS = armaMatE - (-5e-9)*armaMat;

    double d = det(armaMatS);

    arma::vec resArma = inv(armaMat)*armaVec;
    std::vector<double> dJ(dim);


            doIterartion(deviceStatus, plotMutex, simulationData, outputData, 0, iteration);

            cathodeFields.clear();
            for (int i = 0; i < nFlows; i++)
            {
                    if (deviceStatus.GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
                    {
                            cathodeFieldsTmp =
    emissionCurrentSolver.CalculateCathodeFields(deviceStatus.GetFlow(i)->GetEmitterDevice(), &particleGridInterface,
    deviceStatus.GetGridData(), particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY(), i, 0,
    deviceStatus.GetFlow(i)->GetMass(), deviceStatus.GetFlow(i)->GetCharge(), p2 - p1); for (int k = 0; k <
    cathodeFieldsTmp.size(); k++)
                            {
                                    cathodeFields.push_back(cathodeFieldsTmp[k]);
                            }
                    }
            }

            for (int i = 0; i < dim; i++)
            {
                    armaVec(i) = (cathodeFields[i] - Emin);

            }

    }*/

    double max = std::abs(armaMat).max();

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            //		if (std::abs(armaMat(i, j)) < 0.05*max && 1!=j)
            //	armaMat(i, j) = 0;
        }
    }

    double              cond = det(armaMat) * det(inv(armaMat));
    std::vector<double> dJ(dim);

    float stepN;

    int iter = 0;

    std::vector<double> currents1(nFlows);
    std::vector<double> currents2(nFlows);

    while (1)
    {

        arma::mat resArma = inv(armaMat) * armaVec;

        stepN = std::min((0.7 + pow(float(iter) / 50.0, 4)), 0.9);

        for (int flowsIter = 0; flowsIter < nFlows; flowsIter++)
        {
            int i = flowsSol[flowsIter];
            for (int j = 0; j < p2; j++)
                dJ[j]  = stepN * resArma(j + p2 * flowsIter);

            currents1[flowsIter] =
                deviceStatus.GetFlow(i)->GetEmitterDevice()->GetParticleSource()->GetEmissionCurrent(1);
            emissionCurrentSolver.ChangePolinom(deviceStatus.GetFlow(i)->GetEmitterDevice(), dJ, p1, p2);
        };

        doIterartion(deviceStatus, plotMutex, simulationData, outputData, 0, iteration);

        cathodeFields.clear();
        float delta = -1;
        for (int flowsIter = 0; flowsIter < nFlows; flowsIter++)
        {
            int i = flowsSol[flowsIter];

            currents2[flowsIter] =
                deviceStatus.GetFlow(i)->GetEmitterDevice()->GetParticleSource()->GetEmissionCurrent(1);

            if (std::abs(currents2[flowsIter] - currents1[flowsIter]) / currents2[flowsIter] > delta)
                delta = std::abs(currents2[flowsIter] - currents1[flowsIter]) / currents2[flowsIter];

            if (deviceStatus.GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
            {
                cathodeFieldsTmp1 = emissionCurrentSolver.CalculateCathodeFields(
                    deviceStatus.GetFlow(i)->GetEmitterDevice(), &particleGridInterface, deviceStatus.GetGridData(),
                    particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY(), i, 0, deviceStatus.GetFlow(i)->GetMass(),
                    deviceStatus.GetFlow(i)->GetCharge(), p2 - p1);
                for (int k = 0; k < cathodeFieldsTmp1.size(); k++)
                {
                    cathodeFields.push_back(cathodeFieldsTmp1[k]);
                }
            }
        }

        if (delta < 0.01)
            break;

        iter++;
        for (int i     = 0; i < dim; i++)
            armaVec(i) = (cathodeFields[i] - Emin) * (cathodeFields[i] - Emin);
    }
    progress = 1;
    /*int res;
    int dimm = dim;

    dgetrf(&dimm, &dimm, A, &dimm, &dimm, &res);
    char sy = 'N';
    int i1 = 1;

    std::vector<double> btmp = b;

    dgetrs(&sy, &dimm, &i1, A, &dimm, &dimm, &b[0], &i1, &res);

    for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
    {
            deviceStatus.GetFlow(i)->GetDynamicsData()->Reset();
    };

    std::vector<double> test(dim);
    test.resize(dim);
    for (int i = 0; i < dim; i++)
    {
            test[i] = 0;
            for (int j = 0; j < dim; j++)
            {
                    test[i] = test[i] + cathodeFieldsNew[j][i] * btmp[j];
            };
    };*/

    //	for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
    //		delete[] W[i];
};

template void SolverPTI2daxsdouble::doIterartion<device2daxsdouble>(device2daxsdouble& deviceStatus,
                                                                    std::mutex&        plotMutex,
                                                                    SimulationData&    simulationData,
                                                                    std::vector<std::vector<DynamicsData*>>& outputData,
                                                                    int flagSaveData, int& iteration);

template void SolverPTI2daxsfloat::doIterartion<device2daxsfloat>(device2daxsfloat& deviceStatus, std::mutex& plotMutex,
                                                                  SimulationData& simulationData,
                                                                  std::vector<std::vector<DynamicsData*>>& outputData,
                                                                  int flagSaveData, int& iteration);

template void SolverPTI2ddouble::doIterartion<device2ddouble>(device2ddouble& deviceStatus, std::mutex& plotMutex,
                                                              SimulationData&                          simulationData,
                                                              std::vector<std::vector<DynamicsData*>>& outputData,
                                                              int flagSaveData, int& iteration);

template void SolverPTI2dfloat::doIterartion<device2dfloat>(device2dfloat& deviceStatus, std::mutex& plotMutex,
                                                            SimulationData&                          simulationData,
                                                            std::vector<std::vector<DynamicsData*>>& outputData,
                                                            int flagSaveData, int& iteration);

template <class ParticleGridInterfaceType, class EmissionCurrentSolverType, class FieldSolver, class PointType>
template <class deviceType>
void SolverPTI<ParticleGridInterfaceType, EmissionCurrentSolverType, FieldSolver, PointType>::doIterartion(
    deviceType& deviceStatus, std::mutex& plotMutex, SimulationData& simulationData,
    std::vector<std::vector<DynamicsData*>>& outputData, int flagSaveData, int& iteration)
{
    int    step  = 0;
    double omega = 1;
    omega        = 0.2;

    std::vector<DynamicsData*> outputDataTmp;

    if (flagSaveData == 1)
    {
        for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
            outputDataTmp.push_back(new DynamicsData());
    }

    int flagSaveDataLoc = 0;

    fieldSolver.FieldSimulate(deviceStatus.gridData, deviceStatus.mesh);
    std::vector<PointType> vTmp = deviceStatus.gridData.V;
    std::vector<PointType> rhoIter(vTmp.size());
    std::vector<PointType> rhoTmp(vTmp.size());

    while (1)
    {
        for (int kk     = 0; kk < rhoIter.size(); kk++)
            rhoIter[kk] = 0;

        for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
        {
            EmptyPlaces[i].clear();
            EmptyPlaces1[i].clear();

            deviceStatus.GetFlow(i)->GetDynamicsData()->clear();

            std::vector<unsigned int> v;
            deviceStatus.GetFlow(i)->GenerateParticles(v, 1, 1, 0);

            particleGridInterface.initBoundaries(deviceStatus.GetFlow(i)->GetboundaryConditions());

            saveIndexes = deviceStatus.GetFlow(i)->GetNewIndexes(v);

            if (flagSaveData == 1 && flagSaveDataLoc == 1)
                outputDataTmp[i]->Init(saveIndexes, deviceStatus.GetFlow(i)->GetDynamicsData()->SpaseSize(),
                                       deviceStatus.GetFlow(i)->GetMass(), deviceStatus.GetFlow(i)->GetCharge(),
                                       sizeof(PointType));

            step = 0;
            while (deviceStatus.GetFlow(i)->GetDynamicsData()->NParticles != 0)
            {
                if (flagSaveData == 1 && flagSaveDataLoc == 1)
                    outputDataTmp[i]->SetData(deviceStatus.GetFlow(i)->GetDynamicsData()->GetData(),
                                              deviceStatus.GetFlow(i)->GetDynamicsData()->Time /
                                                  (commtools::LIGHT_VELOCITY()));

                EmptyPlaces1[i] = particleGridInterface.InCell(deviceStatus.GetFlow(i)->GetDynamicsData());

                for (int kk = 0; kk < EmptyPlaces1[i].size(); kk++)
                {
                    deviceStatus.GetFlow(i)->GetDynamicsData()->cellsNumbers[EmptyPlaces1[i][kk]] =
                        deviceStatus.GetFlow(i)->GetDynamicsDataTmp()->cellsNumbers[EmptyPlaces1[i][kk]];
                };

                for (int k1 = 0; k1 < deviceStatus.GetFlow(i)->GetDynamicsData()->NParticles; k1++)
                {
                    for (int k2 = 0; k2 < 9; k2++)
                    {
                        Wtmp[i][k1][k2] = W[i][k1][k2];
                    }
                }

                particleGridInterface.Wcalculate(deviceStatus.GetFlow(i)->GetDynamicsData(), W[i]);

                plotMutex.lock();
                if (step != 0)
                    particleGridInterface.Particles2GridPTI(deviceStatus.GetFlow(i)->GetDynamicsDataTmp(), Wtmp[i],
                                                            deviceStatus.GetFlow(i)->GetDynamicsData(), W[i],
                                                            &rhoIter[0],
                                                            particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY());
                plotMutex.unlock();

                deviceStatus.GetFlow(i)->GetDynamicsData()->removeParticle(EmptyPlaces[i]);

                if (flagSaveData == 1 && flagSaveDataLoc == 1)
                    outputDataTmp[i]->SetRemovePTI(EmptyPlaces[i]);

                if (deviceStatus.GetFlow(i)->GetDynamicsData()->NParticles <= 0)
                    break;

                particleGridInterface.Grid2Particles(deviceStatus.GetFlow(i)->GetDynamicsData(),
                                                     deviceStatus.GetGridData(), W[i]);

                particleMover.updateMomentums(deviceStatus.GetFlow(i)->GetDynamicsData(),
                                              deviceStatus.GetFlow(i)->GetAlpha(), i);

                deviceStatus.GetFlow(i)->CopyDynamicsDataToTmp();

                particleMover.updatePositions(deviceStatus.GetFlow(i)->GetDynamicsData(), i);

                plotMutex.lock();
                EmptyPlaces[i] = particleGridInterface.CheckParticlesBoundaries(
                    deviceStatus.boundaries, deviceStatus.conductorList, deviceStatus.GetFlow(i)->GetDynamicsDataTmp(),
                    deviceStatus.GetFlow(i)->GetDynamicsData(), particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY(),
                    deviceStatus.GetFlow(i)->GetCharge(), deviceStatus.GetFlow(i)->GetMass());

                plotMutex.unlock();

                step++;
            }
        }

        particleGridInterface.Charge2Density(&rhoIter[0]);

        if (iteration == 1)
            //	omega = 0.9;

            rhoTmp = deviceStatus.gridData.rho;

        for (int k                       = 0; k < deviceStatus.gridData.rho.size(); k++)
            deviceStatus.gridData.rho[k] = (1 - omega) * deviceStatus.gridData.rho[k] + omega * rhoIter[k];

        fieldSolver.FieldSimulate(deviceStatus.gridData, deviceStatus.mesh);

        if (emissionCurrentSolver.GetParameters() == 0)
        {
            for (int i = 0; i < deviceStatus.GetNumberParticlesFlows(); i++)
            {
                emissionCurrentSolver.UpdateEmissionCurrent(
                    deviceStatus.GetFlow(i)->GetEmitterDevice(), &particleGridInterface, deviceStatus.GetGridData(),
                    particleMover.GetTimeStep(i) / commtools::LIGHT_VELOCITY(), i, step, deviceStatus.GetFlow(i)->GetMass(),
                    deviceStatus.GetFlow(i)->GetCharge());
            }
        }

        PointType iterationError = -1;
        for (int k = 0; k < deviceStatus.gridData.rho.size(); k++)
        {
            vTmp[k] = std::abs((vTmp[k] - deviceStatus.gridData.V[k]) / deviceStatus.gridData.V[k]);
            if (vTmp[k] > iterationError)
                iterationError = vTmp[k];
        }

        vTmp = deviceStatus.gridData.V;

        if (iteration > 0)
        {
            if (simulationData.YData[0][iteration - 1] < iterationError)
            {
                omega = omega * 0.8;
            };
        }

        iteration++;

        simulationData.addData(deviceStatus, iteration, iterationError);

        if (iterationError < 0.00005)
        {
            if (flagSaveDataLoc == 0)
                flagSaveDataLoc = 1;
            else
                break;
        }
    };
    if (flagSaveData == 1)
        outputData.push_back(outputDataTmp);
}
