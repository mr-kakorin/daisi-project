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
#include "Model2dpolar.h"
#include "Particle.h"
#include "ParticleGridInterface.h"
#include "ParticlesFlow.h"
#include "ParticlesMover.h"
#include "PoissonSolver.h"
#include "Results.h"
#include "armadillo"

#include <iostream>
// #include <windows.h>

template Solver<float>;
template Solver<double>;

std::vector<double> testV;

void vectorSubtraction(std::vector<unsigned int>& v1, std::vector<unsigned int> v2)
{
    for (int i    = 0; i < v2.size(); i++)
        v1[v2[i]] = -1;

    int k = int(v1.size());
    for (int i = 0; i < k; i++)
    {
        if (v1[i] == -1)
        {
            v1.erase(v1.begin() + i);
            i--;
            k--;
        };
    };
};

template void Solver<double>::SimulateCPUPTI<device2daxsdouble>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2daxsdouble>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<float>::SimulateCPUPTI<device2daxsfloat>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2daxsfloat>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<double>::SimulateCPUPTI<device2ddouble>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2ddouble>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<float>::SimulateCPUPTI<device2dfloat>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2dfloat>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<double>::SimulateCPUPTI<device2dpolardouble>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2dpolardouble>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<float>::SimulateCPUPTI<device2dpolarfloat>(
    double& progress, bool& flagAbort, const std::shared_ptr<device2dpolarfloat>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<double>::SimulateCPUPTI<device3dExtrdouble>(
    double& progress, bool& flagAbort, const std::shared_ptr<device3dExtrdouble>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template void Solver<float>::SimulateCPUPTI<device3dExtrfloat>(
    double& progress, bool& flagAbort, const std::shared_ptr<device3dExtrfloat>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

template <class PointType>
template <class deviceType>
void Solver<PointType>::SimulateCPUPTI(
    double& progress, bool& flagAbort, const std::shared_ptr<deviceType>& deviceStatus,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData)
{
    double progressLoc;
    double timeDilatation   = parameters[0];
    double outTime          = parameters[1];
    int    saveParam        = parameters[2];
    int    numberSaveTraces = parameters[3];
    int    numThreads       = parameters[4];
    int    blockSize        = int(parameters[6]);

    std::vector<int> nonZeros;

    omp_set_num_threads(numThreads);

    double dH = deviceStatus->GetGridData()->GetMaxSixe();

    outTime = outTime * (1e-9) * commtools::LIGHT_VELOCITY();

    progress = 0.015;

    int memorySize = 0;
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        if (deviceStatus->GetFlow(i)->GetMaxParticles() > memorySize)
            memorySize = deviceStatus->GetFlow(i)->GetMaxParticles();
    }
    deviceStatus->GetGridData()->InitParallel(numThreads);
    particleMover->InitParallel(numThreads);
    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
        deviceStatus->GetElectrodes()[i]->InitParallel(numThreads);

    startCellNumbersGridTmp.resize(deviceStatus->GetNumberParticlesFlows());
    startCellNumbersGridTmp1.resize(deviceStatus->GetNumberParticlesFlows());

    startCellNumbersTmp.resize(deviceStatus->GetNumberParticlesFlows());
    startCellNumbersTmp1.resize(deviceStatus->GetNumberParticlesFlows());

    if (parameters[7] == 1 || fieldSolver->solverFlags[2] == 1)
        InitFieldSolver(deviceStatus, progressLoc);

    if (parameters[7] == 1)
    {
        fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                                   deviceStatus->GetboundaryConditions(),
                                   deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                                       (commtools::LIGHT_VELOCITY()),
                                   progressLoc);
    }

    simulationData.reset();

    std::vector<int> boundaryPoints = fieldSolver->getNearBoundarypoints_all();

    particleGridInterface->init(
        deviceStatus->GetNumberParticlesFlows(), deviceStatus->GetGridData(),
        deviceStatus->Getmesh()->templNumb, deviceStatus->Getmesh()->flagMatrix, boundaryPoints,
        deviceStatus->GetDomainBoundary(), memorySize, numThreads, blockSize);

    EmptyPlaces.resize(numThreads);
    EmptyPlaces1.resize(deviceStatus->GetNumberParticlesFlows());

    std::vector<int> flag;

    emissionCurrentSolverPTI->reset();

    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
    {
        deviceStatus->GetElectrodes()[i]->ResetPower();
    }
    int iteration = 0;

    std::vector<PointType> alpas(deviceStatus->GetNumberParticlesFlows());

    std::vector<int> emTypes(deviceStatus->GetNumberParticlesFlows());

    if (parameters[7] == 1)
    {
        fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                                   deviceStatus->GetboundaryConditions(),
                                   deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                                       (commtools::LIGHT_VELOCITY()),
                                   progressLoc);
    }

    std::vector<unsigned int> v;

    PointType Lmax = 0;
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        PointType l = deviceStatus->GetFlow(i)->GetEmitterDevice()->GetSourceSize();
        if (l > Lmax)
            Lmax = l;
    }

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        emTypes[i] = deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType();
        deviceStatus->GetFlow(i)->ReserveMemory(memorySize);
        emissionCurrentSolverPTI->init(
            i, deviceStatus->GetGridData(), emTypes[i],
            deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSource(),
            particleGridInterface, deviceStatus->GetGridData()->GetType());
        deviceStatus->GetFlow(i)->InitParallel(memorySize, numThreads, 0);
    }

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        flag.push_back(1);
        startCellNumbersGridTmp[i].resize(numThreads);
        startCellNumbersGridTmp1[i].resize(numThreads);

        startCellNumbersTmp[i].resize(numThreads);
        startCellNumbersTmp1[i].resize(numThreads);

        for (int thread = 0; thread < numThreads; thread++)
        {

            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->clear();
            deviceStatus->GetFlow(i)->GenerateParticlesThreaded(thread, numThreads, v, 0, 1, 0,
                                                                deviceStatus->GetGridData(), 0, 1);

            particleGridInterface->SearchStartCells(
                i, deviceStatus->GetFlow(i)->GetDynamicsData(thread));
            startCellNumbersTmp1[i][thread] =
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->cellsNumbers;
            particleGridInterface->InitEmCells(
                emissionCurrentSolverPTI->nearCathodeVolumes,
                deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(),
                deviceStatus->GetFlow(i)->GetDynamicsData(thread), dH,
                deviceStatus->GetGridData()->GetType(), 0);
            startCellNumbersGridTmp1[i][thread] =
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Get_startCellNumbers();

            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->clear();
            deviceStatus->GetFlow(i)->GenerateParticlesThreaded(thread, numThreads, v, 0, 1, 0,
                                                                deviceStatus->GetGridData(), 0, 0);

            particleGridInterface->SearchStartCells(
                i, deviceStatus->GetFlow(i)->GetDynamicsData(thread));
            startCellNumbersTmp[i][thread] =
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->cellsNumbers;
            particleGridInterface->InitEmCells(
                emissionCurrentSolverPTI->nearCathodeVolumes,
                deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(),
                deviceStatus->GetFlow(i)->GetDynamicsData(thread), dH,
                deviceStatus->GetGridData()->GetType(), 1);
            startCellNumbersGridTmp[i][thread] =
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Get_startCellNumbers();
        }

        alpas[i] = deviceStatus->GetFlow(i)->GetAlpha();
    }

    //	int minBlockSize = 2000;
    //	if (blockSize < minBlockSize)
    //		blockSize = minBlockSize;

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
        deviceStatus->GetFlow(i)->InitParallelTmpBlocks(numThreads, blockSize, 1);

    W.resize(numThreads);
    for (int i = 0; i < numThreads; i++)
        W[i]   = new PointType[blockSize][9];

    Wtmp.resize(numThreads);
    for (int i  = 0; i < numThreads; i++)
        Wtmp[i] = new PointType[blockSize][9];

    particleMover->init(alpas, 1);

    //	std::vector<unsigned int> emissionCells = emissionCurrentSolverPTI->GetEmissionCells();

    float maxEr;
    //	errors = estimateError(deviceStatus, maxEr);
    particleMover->flagInit = 0;

    std::vector<double> dt;

    TimeStepEstimate(dt, deviceStatus, numThreads, progressLoc);

    particleMover->flagInit = 1;
    particleMover->SetTimeSteps(dt);

    // doIterartion(flagAbort, emissionCells, 1.0e-4, deviceStatus, plotMutex, simulationData,
    // outputData, 1, iteration,
    // saveParam, numberSaveTraces, numThreads, blockSize);

    progress = 1;
    return;
};

template void Solver<double>::doIterartion<device2daxsdouble>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2daxsdouble>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<float>::doIterartion<device2daxsfloat>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2daxsfloat>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<double>::doIterartion<device2ddouble>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2ddouble>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<float>::doIterartion<device2dfloat>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2dfloat>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<double>::doIterartion<device2dpolardouble>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2dpolardouble>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<float>::doIterartion<device2dpolarfloat>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device2dpolarfloat>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<double>::doIterartion<device3dExtrdouble>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device3dExtrdouble>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template void Solver<float>::doIterartion<device3dExtrfloat>(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<device3dExtrfloat>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize);

template <class PointType>
template <class deviceType>
void Solver<PointType>::doIterartion(
    bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
    const std::shared_ptr<deviceType>& deviceStatus, std::mutex& plotMutex,
    const std::shared_ptr<SimulationData>&                                simulationData,
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
    int flagSaveData, int& iteration, int saveParam, int numberSaveTraces, int numThreads,
    int blockSize)
{

    std::vector<double> ChargeSign;
    double              progressLoc;
    double              omega  = relaxations[0];
    double              omegaB = relaxations[1];

    int                        flagStopUpdateCurrent = 0;
    std::vector<DynamicsData*> outputDataTmp;

    if (flagSaveData == 1)
    {
        for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
            outputDataTmp.push_back(new DynamicsData());
    }
    std::vector<int> nonZeros;
    nonZeros.reserve(1000);

    int flagSaveDataLoc = 0;

    // fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
    // deviceStatus->GetboundaryConditions(), 0);
    std::vector<PointType> vTmp = deviceStatus->GetGridData()->GetV();
    std::vector<PointType> rhoTmp(vTmp.size());
    std::vector<PointType> BTmp;  //= deviceStatus->GetGridData()->GetB();
    std::vector<PointType> BTmp1; //= deviceStatus->GetGridData()->GetB();

    std::vector<PointType> currentsTmp;

    int flagDistrib = 1;

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        ChargeSign.push_back(deviceStatus->GetFlow(i)->GetCharge());

        if (deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
            currentsTmp.push_back(
                deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionCurrent());
    }

    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
    {
        deviceStatus->GetElectrodes()[i]->ResetPower();
        deviceStatus->GetElectrodes()[i]->ResetCurrent();
    };

    /*flagDistrib = 1;
    flagStopUpdateCurrent = 1;
    flagAbort = false;
    flagSaveDataLoc = 1;*/

    int totalParticels;

    std::vector<arma::mat> Icoef(numThreads);
    volatile int           flagAllBreak = 0;

/*std::string currentModelFilePath = "C:/gesaSlice1/modelFiles/gesaSlice1_0.mdl";
std::ifstream ifs(currentModelFilePath.c_str());
boost::archive::binary_iarchive ia(ifs);

Model2dpolarPTIdouble* ser2 = new Model2dpolarPTIdouble();
ia >> *ser2;

std::vector<int> bp;

ser2->solver.InitfieldSolver(ser2->deviceStatus);
//ser2->solver.fieldSolver->FieldSimulate(ser2->deviceStatus->GetGridData(),
ser2->deviceStatus->Getmesh(),
*ser2->deviceStatus->GetboundaryConditions(), 0);


ser2->solver.particleGridInterface->init(ser2->deviceStatus->GetGridData(),
ser2->deviceStatus->Getmesh()->templNumb,
ser2->deviceStatus->Getmesh()->flagMatrix, bp, ser2->deviceStatus->GetDomainBoundary(), 1e5,
numThreads);
*/

#pragma omp parallel
    {

        while (1)
        {
#pragma omp barrier
#pragma omp single
            {

                for (int thread = 1; thread < numThreads; thread++)
                    Icoef[0]    = Icoef[0] + Icoef[thread];

                nonZeros.clear();
                deviceStatus->GetGridData()->Summrho();
                particleGridInterface->Charge2Density(&deviceStatus->GetGridData()->Getrho()[0],
                                                      nonZeros);

                magneticFieldSolver->FieldSimulate(
                    deviceStatus->GetElectrodes(), deviceStatus->GetGridData(),
                    particleMover->GetTimeStep(0, 0) / commtools::LIGHT_VELOCITY());

                //	BTmp1 = deviceStatus->GetGridData()->GetB();

                for (int k = 0; k < deviceStatus->GetGridData()->Getrho().size(); k++)
                {
                    deviceStatus->GetGridData()->Getrho()[k] =
                        (1 - omega) * rhoTmp[k] + omega * deviceStatus->GetGridData()->Getrho()[k];
                }

                for (int k = 0; k < BTmp.size(); k++)
                {
                    BTmp[k] = (1 - omegaB) * BTmp[k] + omegaB * BTmp1[k];
                }

                //		deviceStatus->GetGridData()->SetB(BTmp);

                //	BTmp = deviceStatus->GetGridData()->GetB();

                // fieldSolver->FieldSimulateCharge(deviceStatus->GetGridData(),
                // deviceStatus->Getmesh(),
                // deviceStatus->GetboundaryConditions(), 0, nonZeros, 0);
                fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                                           deviceStatus->GetboundaryConditions(), 0, progressLoc);
                // deviceStatus->GetGridData()->ApplyTimeDepending(deviceStatus->LinacParams[0], 0,
                // deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                // (commtools::LIGHT_VELOCITY()), nonZeros);

                for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
                {
                    deviceStatus->GetElectrodes()[i]->ResetCurrent();
                    //	deviceStatus->GetElectrodes()[i]->PowerAndCurrentsCalculate(0,1);
                }
                int flagE = 0;

                if (emissionCurrentSolverPTI->GetParameters()[0][0] < 2 && iteration > 0)
                {
                    if (emissionCurrentSolverPTI->GetParameters()[0][0] == 0)
                    {
                        //			emissionCurrentSolverPTI->UpdateEmissionCurrent(flagE,
                        //Icoef[0],
                        // deviceStatus->GetEmittersVector(), &particleGridInterface,
                        // deviceStatus->GetGridData(),
                        // flagStopUpdateCurrent, ChargeSign);
                    }
                }

                int       ii = 0;
                PointType currentNew;
                PointType errorCurrent = -1;
                for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
                {
                    if (emissionCurrentSolverPTI->GetParameters()[0][0] < 2 && iteration > 0)
                    {
                        if (iteration % 1 == 0 && flagStopUpdateCurrent == 0 &&
                            emissionCurrentSolverPTI->GetParameters()[0][0] == 1)
                        {
                            //	emissionCurrentSolverPTI->UpdateEmissionCurrent(flagE,
                            // deviceStatus->GetFlow(i)->GetEmitterDevice(), &particleGridInterface,
                            // deviceStatus->GetGridData(), particleMover->GetTimeStep(i, 0) /
                            // commtools::LIGHT_VELOCITY(),
                            // i,  1, deviceStatus->GetFlow(i)->GetMass(),
                            // deviceStatus->GetFlow(i)->GetCharge());
                        }
                    }

                    if (deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
                    {
                        currentNew =
                            deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionCurrent();
                        PointType err = std::abs((currentsTmp[ii] - currentNew) / currentNew);
                        if (errorCurrent < err)
                            errorCurrent = err;

                        currentsTmp[ii] = currentNew;
                        ii++;
                    }
                }

                std::vector<unsigned int> v;
                std::vector<unsigned int> vDel;

                for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
                {
                    for (int thread = 0; thread < numThreads; thread++)
                    {
                        vDel.clear();
                        deviceStatus->GetFlow(i)->GenerateParticlesThreaded(
                            thread, numThreads, v, 0, 1, 0, deviceStatus->GetGridData(), 0,
                            flagDistrib);

                        if (flagDistrib == 0)
                        {
                            deviceStatus->GetFlow(i)
                                ->GetDynamicsData(thread)
                                ->Get_startCellNumbers() = startCellNumbersGridTmp[i][thread];
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->cellsNumbers =
                                startCellNumbersTmp[i][thread];
                            //		deviceStatus->GetFlow(i)->GetDynamicsData(thread)->startCellNumbers
                            //=
                            // startCellNumbersTmp[i][thread];
                        }
                        else
                        {
                            deviceStatus->GetFlow(i)
                                ->GetDynamicsData(thread)
                                ->Get_startCellNumbers() = startCellNumbersGridTmp1[i][thread];
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->cellsNumbers =
                                startCellNumbersTmp1[i][thread];
                            //	deviceStatus->GetFlow(i)->GetDynamicsData(thread)->startCellNumbers
                            //=
                            // startCellNumbersTmp[i][thread];
                        }
                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->removeParticle(vDel);
                    }
                }

                if (errorCurrent < 0.0001 && iteration > 2)
                {
                    flagStopUpdateCurrent = 1;
                    omega                 = relaxations[2];
                    omegaB                = relaxations[3];
                };

                PointType iterationError = -1;
                for (int k = 0; k < deviceStatus->GetGridData()->Getrho().size(); k++)
                {
                    if (std::abs(deviceStatus->GetGridData()->GetV()[k]) < 0.01)
                        vTmp[k] = 0;
                    else
                        vTmp[k] = std::abs((vTmp[k] - deviceStatus->GetGridData()->GetV()[k]) /
                                      deviceStatus->GetGridData()->GetV()[k]);
                    if (vTmp[k] > iterationError)
                        iterationError = vTmp[k];
                }

                vTmp = deviceStatus->GetGridData()->GetV();

                /*if (iteration > 0)
                {
                if (simulationData.YData[0][iteration - 1] < iterationError && flagStopUpdateCurrent
                == 1)
                {
                omega = omega*0.7;
                };
                }*/

                iteration++;

                //	simulationData.addData(deviceStatus, iteration, iterationError);

                /*	if (flagE)
                {
                deviceStatus->GetGridData()->Getrho() = rhoTmp;
                fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(), 0);
                };*/

                //	iterationError = 1e-5;
                if (iterationError < tolerance && 0 == flagDistrib && iteration > 5)
                {
                    flagDistrib           = 1;
                    flagStopUpdateCurrent = 1;
                    flagAbort             = false;
                    flagSaveDataLoc       = 1;
                }
                else
                {
                    volatile int flagBreak = 0;
                    if ((iterationError < tolerance && iteration > 5) || flagAbort == false)
                    {
                        if (0 == flagSaveDataLoc)
                        {
                            flagSaveDataLoc       = 1;
                            flagDistrib           = 1;
                            flagStopUpdateCurrent = 1;
                            flagBreak             = 1;
                        }
                    }
                    if (1 == flagSaveDataLoc && 0 == flagBreak)
                        flagAllBreak = 1;
                }

                for (int thread = 0; thread < numThreads; thread++)
                    Icoef[thread].zeros(emissionCurrentSolverPTI->GetEmSize(),
                                        emissionCurrentSolverPTI->GetEmSize());

                rhoTmp = deviceStatus->GetGridData()->Getrho();
                deviceStatus->GetGridData()->densityReset();

                for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
                {
                    std::vector<unsigned int> v;

                    emissionCurrentSolverPTI->CalculateCathodeFields(
                        deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSource(),
                        deviceStatus->GetGridData(), i);

                    totalParticels = deviceStatus->GetFlow(i)->GetNumberOfParticles();

                    if (flagSaveData == 1 && flagSaveDataLoc == 1)
                        outputDataTmp[i]->Init(
                            deviceStatus->GetFlow(i)->GetDynamicsData()->SpaseSize(),
                            deviceStatus->GetFlow(i)->GetMass(),
                            deviceStatus->GetFlow(i)->GetCharge(), sizeof(PointType),
                            std::min(totalParticels, numberSaveTraces), numThreads, 0);

                    deviceStatus->GetFlow(i)->CalculateFlowCurrent();

                    if (deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType() == 0)
                        particleGridInterface->SearchStartCellsEmission(
                            emissionCurrentSolverPTI->nearCathodeVolumes,
                            deviceStatus->GetFlow(i)->GetDynamicsDataParallelArray());
                }
            }

            if (flagAllBreak == 1)
                break;

            std::vector<unsigned int> v;

            for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
            {
                int thread = omp_get_thread_num();

                EmptyPlaces[thread].clear();
                int i1;
                int i2;
                int step           = 0;
                int blockSizeLocal = std::min(
                    blockSize, deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles());

                int nBlocks;
                if (blockSizeLocal == 0)
                    nBlocks = 0;
                else
                    nBlocks = deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles() /
                              blockSizeLocal;

                std::vector<int> blockIndexes;

                blockIndexes.push_back(0);
                for (int i = 0; i < nBlocks - 1; i++)
                    blockIndexes.push_back(blockSizeLocal * (i + 1));

                blockIndexes.push_back(
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles());

                for (int p = 1; p < blockIndexes.size(); p++)
                {

                    i1 = 0;
                    i2 = blockIndexes[p] - blockIndexes[p - 1];

                    if (i1 == i2)
                        break;

                    step = 0;

                    std::vector<unsigned int> indexes;
                    for (int s = i1; s < i2; s++)
                        indexes.push_back(s);

                    if (flagSaveData == 1 && flagSaveDataLoc == 1)
                        outputDataTmp[i]->AddBlock(indexes, thread, nBlocks, p - 1);

                    deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time = 0;

                    EmptyPlaces[thread].clear();

                    indexes.reserve(i2 - i1);
                    while (i2 != 0)
                    {
                        /*indexes.clear();
                        indexes.resize(i2-i1);
                        for (int s = i1; s < i2; s++)
                        indexes.push_back(s);*/

                        if (flagSaveData == 1 && flagSaveDataLoc == 1 && 0 == step % saveParam)
                            outputDataTmp[i]->SetData(
                                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->GetData(),
                                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time /
                                    (commtools::LIGHT_VELOCITY()),
                                thread, 1, 0);

                        if (step < 2)
                            particleGridInterface->InCellWithEps(
                                deviceStatus->GetFlow(i)->GetDynamicsData(thread), i1, i2, thread,
                                i);
                        else
                            particleGridInterface->InCell(
                                deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                                EmptyPlaces[thread], i1, i2, thread, i);

                        //	ser2->solver.particleGridInterface->axsPolar(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                        // ser2->deviceStatus->GetGridData(), EmptyPlaces[thread], i1, i2, thread);

                        /*	for (int k2 = 0; k2 < 4; k2++)
                                {
                                for (int k1 = 0; k1 < indexes.size(); k1++)
                                Wtmp[thread][k1][k2] = W[thread][indexes[k1]][k2];
                                }*/

                        particleGridInterface->Wcalculate(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), W[thread], i1, i2,
                            thread, particleMover->GetParameters()[0]);

                        if (step > 1)
                            particleGridInterface->Particles2GridPTI(
                                emissionCurrentSolverPTI->nearCathodeVolumes,
                                deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(),
                                Icoef[thread], emissionCells,
                                deviceStatus->GetFlow(i)->GetDynamicsDataTmp(thread), Wtmp[thread],
                                deviceStatus->GetFlow(i)->GetDynamicsData(thread), W[thread],
                                &deviceStatus->GetGridData()->Getrho(thread)[0],
                                particleMover->GetTimeStep(i, thread) / commtools::LIGHT_VELOCITY(),
                                deviceStatus->GetGridData()->GetType(), i1, i2, thread);
                        // particleGridInterface->Particles2GridPTI(emissionCurrentSolverPTI->nearCathodeVolumes,
                        // deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(), Icoef,
                        // emissionCells,
                        // deviceStatus->GetFlow(i)->GetDynamicsDataTmp(), Wtmp[i],
                        // deviceStatus->GetFlow(i)->GetDynamicsData(), W[i], &rhoIter[0],
                        // particleMover->GetTimeStep(i)
                        // / commtools::LIGHT_VELOCITY(), deviceStatus->GetGridData()->GetType());

                        particleGridInterface->Grid2Particles(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), fields[thread],
                            deviceStatus->GetGridData(), W[thread], i1, i2, step, recalculate);

                        /*if (flagSaveData == 1 && flagSaveDataLoc == 1)
                                outputDataTmp[i]->SetRemovePTI(EmptyPlaces[thread], thread);

                        i2 = i2 - EmptyPlaces[thread].size();

                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->removeParticle(EmptyPlaces[thread]);
                        EmptyPlaces[thread].clear();*/

                        if (step % 4 == 0)
                        {
                            i2 = i2 - EmptyPlaces[thread].size();
                            if (flagSaveData == 1 && flagSaveDataLoc == 1)
                                outputDataTmp[i]->SetRemovePTI(EmptyPlaces[thread], thread);
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->removeParticle(
                                EmptyPlaces[thread]);
                            EmptyPlaces[thread].clear();
                        }

                        //

                        particleGridInterface->Wcalculate(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), Wtmp[thread], i1, i2,
                            thread, 0);

                        volatile int np =
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles();

                        if (i1 == i2 || 0 == np)
                            break;

                        //	vectorSubtraction(indexes, EmptyPlaces[thread]);

                        deviceStatus->GetFlow(i)->CopyDynamicsDataToTmpThreaded(thread, i1, i2);

                        // particleGridInterface->Wcalculate(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                        // Wtmp[thread], i1, i2);

                        particleMover->stepEstimate(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), i, i1, i2, thread,
                            0);

                        //		testV.push_back(deviceStatus->GetFlow(i)->GetDynamicsData(thread)->GetPointerToPosition1()[0]);

                        if (step == 250)
                        {
                            int tt = 0;
                        };

                        particleMover->updateMomentums(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), fields[thread],
                            deviceStatus->GetFlow(i)->GetAlpha(), i, i1, i2, thread);
                        particleMover->updatePositions(
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                            deviceStatus->GetFlow(i)->GetDynamicsDataTmp(thread), i, i1, i2,
                            thread);

                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time =
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time +
                            particleMover->GetTimeStep(i, thread);

                        particleGridInterface->CheckParticlesBoundaries(
                            deviceStatus->GetFlow(i)->GetboundaryConditions(), EmptyPlaces[thread],
                            deviceStatus->GetboundariesForFlows(), deviceStatus->GetElectrodes(),
                            deviceStatus->GetFlow(i)->GetDynamicsDataTmp(thread),
                            deviceStatus->GetFlow(i)->GetDynamicsData(thread), 1,
                            deviceStatus->GetFlow(i)->GetCharge(),
                            deviceStatus->GetFlow(i)->GetMass(), i1, i2, thread);
                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->searchBadParticle(
                            EmptyPlaces[thread]);

                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->markBadParticles(
                            EmptyPlaces[thread]);

                        if (deviceStatus->GetFlow(i)->CheckTimeLimit(thread))
                        {
                            EmptyPlaces[thread].clear();
                            for (int kk = 0;
                                 kk <
                                 deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles();
                                 kk++)
                                EmptyPlaces[thread].push_back(kk);
                        };

                        step++;
                    }
                }
            }
        }
    }
    //	if (flagSaveData == 1)
    //	outputData.push_back(outputDataTmp);
}
