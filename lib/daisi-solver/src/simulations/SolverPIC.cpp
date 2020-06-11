#include "Solver.h"
#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmissionCurrentSolver.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice3d.h"
#include "FlagStringsSolver.h"
#include "GridData.h"
#include "MagneticFieldSolver.h"
#include "MeshGenerator.h"
#include "Particle.h"
#include "ParticleGridInterface.h"
#include "ParticlesFlow.h"
#include "ParticlesMover.h"
#include "PoissonSolver.h"
#include "Results.h"
#include <iostream>
#include <omp.h>


template <class PointType>
template <class deviceType>
void Solver<PointType>::SimulateCPUPIC(
    double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
    const std::shared_ptr<deviceType>&                                    deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
    std::string& errorMsg)
{
    progressLoc = 0;
    status.push_back("PIC simulations process started");
    status.push_back("Checking configuration");

    CheckConfiguration(deviceStatus, flagRestart, errorMsg);
    if (errorMsg.size())
        return;
    status.push_back("Memory initialization and allocation");

    InitLocalVariables(deviceStatus->GetNumberParticlesFlows(),
                       deviceStatus->GetFlow(0)->GetDynamicsData()->getFieldSize());

    InitOutputData(outputData, deviceStatus->GetNumberParticlesFlows(),
                   deviceStatus->GetFlow(0)->GetDynamicsData()->SpaseSize(),
                   deviceStatus->GetMasses(), deviceStatus->GetCharges(), sizeof(PointType));

    deviceStatus->initParallel(numThreads, flagRestart, memorySize, 1, blockSize);

    status.push_back("Solvers initialization");
    InitSolvers(deviceStatus, flagRestart, progressLoc, status);

    int step = -1;

    volatile int flagAllBreak = 0;

    status.push_back("Main PIC loop started");

#pragma omp parallel num_threads(numThreads)
    {

        while (1)
        {

#pragma omp barrier
#pragma omp single
            {
                deviceStatus->SyncronizeThreads(step, particleMover->GetTimeStep(0, 0) /
                                                          (LIGHT_VELOCITY() * (1e-9)));

                plotMutex.lock();
                nonZeros.clear();
                particleGridInterface->Charge2Density(&deviceStatus->GetGridData()->Getrho()[0],
                                                      nonZeros);
                if (nonZeros.size())
                {
                    if (deviceStatus->isTimeVarying())
                        fieldSolver->FieldSimulate(
                            deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                            deviceStatus->GetboundaryConditions(),
                            deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                                (LIGHT_VELOCITY()),
                            progressLoc);

                    fieldSolver->FieldSimulateCharge(
                        deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                        deviceStatus->GetboundaryConditions(),
                        deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                            (LIGHT_VELOCITY()),
                        nonZeros, step);
                    //	magneticFieldSolver->FieldSimulate(deviceStatus->GetElectrodes(),
                    // deviceStatus->GetGridData(),
                    // particleMover->GetTimeStep(0, 0) / LIGHT_VELOCITY());
                    deviceStatus->GetGridData()->ApplyTimeDepending(
                        deviceStatus->GetglobalFieldConditions(),
                        deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                            (LIGHT_VELOCITY()));
                }
                plotMutex.unlock();

                for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
                {
                    for (int thread = 0; thread < numThreads; thread++)
                    {
                        for (int t = 0; t < outTimesLoc.size(); t++)
                        {
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->GetSaveIndexes(
                                saveIndexes, outTimesLoc[t], step, saveParam, tracesSaveProbability,
                                flagSave[i][thread][t]);
                            outputData.back()[t][i]->AddSavedTraces(saveIndexes, thread);
                        }
                    }
                }

                for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
                    emissionCurrentSolverPIC->UpdateEmissionCurrent(
                        deviceStatus->GetFlow(i)->GetEmitterDevice(), particleGridInterface,
                        deviceStatus->GetGridData(),
                        particleMover->GetTimeStep(i, 0) / LIGHT_VELOCITY(), i, step,
                        deviceStatus->GetFlow(i)->GetMass(), deviceStatus->GetFlow(i)->GetCharge(),
                        deviceStatus->GetFlow(i)->GetDistributionStyle());

                if (step > 1 && (step % 2 == 0 || step == 0 || step == 1))
                    simulationData->addDataPIC(deviceStatus);

                step++;
                progress = deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                           (maxtime * 1e-9 * LIGHT_VELOCITY());
                progressLoc = progress;
                //	if (step>1 && (flagAbort == false || (iterationErrorEr<1.0e-9 && step>300)
                //||
                // deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time>particleMover->GetTimeLimit()))

                //volatile double tt = deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time;
                if ((step > 1 && (!flagAbort ||
                                  deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time >
                                      maxtime * 1e-9 * LIGHT_VELOCITY())))
                    flagAllBreak = 1;
                else
                {
                    //	if (fieldSolver->solverFlags[2] == 1)
                    deviceStatus->GetGridData()->densityReset();
                    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
                        deviceStatus->GetElectrodes()[i]->ResetCurrent();
                }
            }

            //#pragma omp barrier

            if (flagAllBreak == 1)
                break;

            for (int iflow = 0; iflow < deviceStatus->GetNumberParticlesFlows(); iflow++)
            {

                volatile int i1;
                volatile int i2      = -1;
                volatile int thread0 = omp_get_thread_num();

                volatile int blockSizeLocal =
                    std::min(blockSize,
                             deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->NParticles());

                volatile int nBlocks;

                if (blockSizeLocal == 0)
                    nBlocks = 1;
                else
                    nBlocks = deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->NParticles() /
                              blockSizeLocal;

                for (volatile int ii          = 0; ii < nBlocks - 1; ii++)
                    blockIndexes[thread0][ii] = blockSizeLocal * (ii + 1);

                blockIndexes[thread0][nBlocks - 1] =
                    deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->NParticles();

                i2 = 0;
                //	if (nonZeros.size())
                //	{

                for (volatile int p = 0; p < nBlocks; p++)
                {

                    i1 = i2;
                    i2 = blockIndexes[thread0][p];

                    // if (i1 == i2)
                    //		break;

                    particleGridInterface->Wcalculate(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0), W[thread0], i1, i2,
                        thread0, 0);

                    particleGridInterface->Grid2Particles(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0), fields[thread0],
                        deviceStatus->GetGridData(), W[thread0], i1, i2, step, recalculate);

                    deviceStatus->GetFlow(iflow)->CopyDynamicsDataToTmpThreaded(thread0, i1, i2);
                    particleMover->stepEstimate(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0), iflow, i1, i2,
                        thread0, 1);

                    particleMover->updateMomentums(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0), fields[thread0],
                        deviceStatus->GetFlow(iflow)->GetAlpha(), iflow, i1, i2, thread0);
                    particleMover->updatePositions(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0),
                        deviceStatus->GetFlow(iflow)->GetDynamicsDataTmp(thread0), iflow, i1, i2,
                        thread0);

                    deviceStatus->GetFlow(iflow)->checkEmittances(thread0, i1, i2);

                    particleGridInterface->CheckParticlesBoundaries(
                        deviceStatus->GetFlow(iflow)->GetboundaryConditions(),
                        EmptyPlacesPIC[iflow][thread0], deviceStatus->GetboundariesForFlows(),
                        deviceStatus->GetElectrodes(),
                        deviceStatus->GetFlow(iflow)->GetDynamicsDataTmp(thread0),
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0),
                        particleMover->GetTimeStep(iflow, 0) / LIGHT_VELOCITY(),
                        deviceStatus->GetFlow(iflow)->GetCharge(),
                        deviceStatus->GetFlow(iflow)->GetMass(), i1, i2, thread0);

                    particleGridInterface->InCell(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0),
                        EmptyPlacesPIC[iflow][thread0], i1, i2, thread0, iflow);

                    particleGridInterface->Wcalculate(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0), W[thread0], i1, i2,
                        thread0, 0);

                    particleGridInterface->Particles2Grid(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0),
                        &deviceStatus->GetGridData()->Getrho(thread0)[0], W[thread0], i1, i2);
                }

                deviceStatus->GetFlow(iflow)->AddEmittancesData(thread0);
                for (int t = 0; t < outTimesLoc.size(); t++)
                {
                    outputData.back()[t][iflow]->SetData(
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->GetData(),
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->Time /
                            (LIGHT_VELOCITY()),
                        thread0, saveParam, 1);
                    outputData.back()[t][iflow]->SetEmptyPlaces(EmptyPlacesPIC[iflow][thread0],
                                                                thread0);
                }
                //}

                if (deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->Time <
                    deviceStatus->GetFlow(iflow)->GetFlowProperties()[4] * 1e-9 *
                        LIGHT_VELOCITY())
                    deviceStatus->GetFlow(iflow)->GenerateParticlesThreaded(
                        thread0, numThreads, EmptyPlacesPIC[iflow][thread0], 1,
                        particleMover->GetTimeStep(iflow, thread0) / LIGHT_VELOCITY(),
                        step, deviceStatus->GetGridData(), 0, 1); // deviceStatus->GetGridData()

                deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->searchBadParticle(
                    EmptyPlacesPIC[iflow][thread0]);

                if (step % 20 == 0)
                {
                    for (int t = 0; t < outTimesLoc.size(); t++)
                        outputData.back()[t][iflow]->SetRemovePTI(EmptyPlacesPIC[iflow][thread0],
                                                                  thread0);
                    deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->removeParticle(
                        EmptyPlacesPIC[iflow][thread0]);
                    EmptyPlacesPIC[iflow][thread0].clear();
                }

                if (step > 0)
                    deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->Time =
                        deviceStatus->GetFlow(iflow)->GetDynamicsData(thread0)->Time +
                        particleMover->GetTimeStep(iflow, thread0);
            }
        }
    }
    progress    = 1;
    progressLoc = 1;

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        deviceStatus->GetFlow(i)->MergeEmittancesData();
        deviceStatus->GetFlow(i)->GetDynamicsData()->Reset();
    }

    for (int i = 0; i < numThreads; i++)
        delete[] W[i];

    restartPossible = 1;
    status.push_back("Done");

    flagAbort = false;
}


template class Solver<float>;
template class Solver<double>;

template void Solver<double>::SimulateCPUPIC<device2daxsdouble>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2daxsdouble>&                             deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<float>::SimulateCPUPIC<device2daxsfloat>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2daxsfloat>&                              deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<double>::SimulateCPUPIC<device2ddouble>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2ddouble>&                                deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<float>::SimulateCPUPIC<device2dfloat>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2dfloat>&                                 deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<double>::SimulateCPUPIC<device2dpolardouble>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2dpolardouble>&                           deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<float>::SimulateCPUPIC<device2dpolarfloat>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device2dpolarfloat>&                            deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<double>::SimulateCPUPIC<device3dExtrdouble>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device3dExtrdouble>&                            deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<float>::SimulateCPUPIC<device3dExtrfloat>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device3dExtrfloat>&                             deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);

template void Solver<double>::SimulateCPUPIC<device3ddouble>(
        double& progress, double& progressLoc, std::vector<std::string>& status, bool& flagAbort,
        const std::shared_ptr<device3ddouble>&                                deviceStatus,
        std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
        std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData, int flagRestart,
        std::string& errorMsg);