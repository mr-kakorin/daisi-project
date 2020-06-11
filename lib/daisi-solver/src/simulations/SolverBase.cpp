#include "Solver.h"
#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmissionCurrentSolver.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice3d.h"
#include "MagneticFieldSolver.h"
#include "MeshGenerator.h"
#include "Particle.h"
#include "ParticleGridInterface.h"
#include "ParticlesFlow.h"
#include "ParticlesMover.h"
#include "PoissonSolver.h"
#include "Results.h"


template <class PointType>
template <class deviceType>
void Solver<PointType>::TimeStepEstimate(std::vector<double>&               result,
                                         const std::shared_ptr<deviceType>& deviceStatus,
                                         int numThreads, double& progress)
{
    /*if (deviceStatus->GetFlow(0)->GetDistributionStyle() == 1)
    {
            result.resize(deviceStatus->GetNumberParticlesFlows());
            for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
            {
                    result[i] = deviceStatus->GetFlow(0)->GetEmitterDevice()->GetLambda() / 100.0;
            }
            return;
    }*/

    std::vector<std::vector<double>> resultTmp;

    resultTmp.resize(numThreads);
    result.resize(deviceStatus->GetNumberParticlesFlows());
    double dH = deviceStatus->GetGridData()->GetMaxSixe();

    int flagSaveDataLoc = 0;

    int step = 0;
    int i1, i2;

    resultTmp.resize(deviceStatus->GetNumberParticlesFlows());
    EmptyPlaces.resize(numThreads);

    std::vector<unsigned int> t;
    double                    allTime = 0;
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        double currentStepDefault;
        if (deviceStatus->GetFlow(i)->GetDistributionStyle() == 5)
            currentStepDefault = deviceStatus->GetFlow(i)->GetEmitterDevice()->GetLambda() / 10;
        else
            currentStepDefault = deviceStatus->GetFlow(i)->GetFlowProperties()[4] * 1e-9 *
                                 LIGHT_VELOCITY() / 100;

        for (int thread = 0; thread < numThreads; thread++)
        {
            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time = 0;

            EmptyPlaces[thread].clear();

            deviceStatus->GetFlow(i)->GetDynamicsData(thread)->clear();
            std::vector<unsigned int> v;

            step               = 0;
            double currentStep = currentStepDefault;
            resultTmp[i].push_back(1e9);
            int flagStep = 0;
            //	while (deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles() != 0)
            while (deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time <
                   deviceStatus->GetFlow(i)->GetFlowProperties()[4] * 1e-9 *
                       LIGHT_VELOCITY())
            {
                deviceStatus->GetFlow(i)->GenerateParticlesThreadedTest(
                    thread, numThreads, v, 0, currentStep / LIGHT_VELOCITY(), step,
                    deviceStatus->GetGridData(), 0, 0);

                progress = allTime / (deviceStatus->GetNumberParticlesFlows() * numThreads *
                                      deviceStatus->GetFlow(i)->GetFlowProperties()[4] * 1e-9 *
                                      LIGHT_VELOCITY());

                //	fieldSolver->FieldSimulate(deviceStatus->GetGridData(),
                // deviceStatus->Getmesh(),
                // deviceStatus->GetboundaryConditions(),
                // deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                //(LIGHT_VELOCITY())); 	if (deviceStatus->LinacParams[0] > 1e-3)
                //		deviceStatus->GetGridData()->ApplyTimeDepending(deviceStatus->LinacParams[0],
                // 0,
                // deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time /
                // (LIGHT_VELOCITY()));
                //	deviceStatus->GetGridData()->ApplyTimeDepending(deviceStatus->LinacParams[0],
                // 0,
                // deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                // (LIGHT_VELOCITY()), nonZeros);

                i1 = 0;
                i2 = deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles();

                if (i1 == i2)
                {
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time =
                        deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time + currentStep;
                    allTime = allTime + currentStep;
                    step++;
                    continue;
                }
                deviceStatus->GetGridData()->ApplyTimeDepending(
                    deviceStatus->GetglobalFieldConditions(),
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time /
                        (LIGHT_VELOCITY()));

                if (step < 2)
                    particleGridInterface->InCellWithEps(
                        deviceStatus->GetFlow(i)->GetDynamicsData(thread), i1, i2, thread, i);
                else
                    particleGridInterface->InCell(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                                                  EmptyPlaces[thread], i1, i2, thread, i);

                particleGridInterface->Wcalculate(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                                                  W[thread], i1, i2, thread, 1);

                //	if (deviceStatus->GetFlow(i)->GetDynamicsData(thread)->NParticles() <= 0)
                //		break;

                particleGridInterface->Grid2Particles(
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread), fields[thread],
                    deviceStatus->GetGridData(), W[thread], i1, i2, step, recalculate);

                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->removeParticle(
                    EmptyPlaces[thread]);
                i2 = i2 - EmptyPlaces[thread].size();

                EmptyPlaces[thread].clear();

                //	if (i1 == i2)
                //		break;

                deviceStatus->GetFlow(i)->CopyDynamicsDataToTmpThreaded(thread, i1, i2);

                particleMover->stepEstimate(deviceStatus->GetFlow(i)->GetDynamicsData(thread), i,
                                            i1, i2, thread, 0);

                particleMover->updatePositions(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                                               deviceStatus->GetFlow(i)->GetDynamicsDataTmp(thread),
                                               i, i1, i2, thread);
                particleMover->updateMomentums(deviceStatus->GetFlow(i)->GetDynamicsData(thread),
                                               fields[thread], deviceStatus->GetFlow(i)->GetAlpha(),
                                               i, i1, i2, thread);

                particleGridInterface->CheckParticlesBoundaries(
                    deviceStatus->GetFlow(i)->GetboundaryConditions(), EmptyPlaces[thread],
                    deviceStatus->GetboundariesForFlows(), deviceStatus->GetElectrodes(),
                    deviceStatus->GetFlow(i)->GetDynamicsDataTmp(thread),
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread), 1,
                    deviceStatus->GetFlow(i)->GetCharge(), deviceStatus->GetFlow(i)->GetMass(), i1,
                    i2, thread);
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->searchBadParticle(
                    EmptyPlaces[thread]);

                currentStep = particleMover->GetTimeStep(i, thread);

                if (currentStep < 0)
                    currentStep = currentStepDefault;

                if (currentStep > 0)
                {
                    flagStep = 1;
                    if (resultTmp[i][thread] > currentStep)
                        resultTmp[i][thread] = currentStep;
                }

                //	deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Get_r()emoveBadParticle();
                deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time =
                    deviceStatus->GetFlow(i)->GetDynamicsData(thread)->Time + currentStep;
                allTime = allTime + currentStep;

                step++;
            }
            if (!flagStep)
                resultTmp[i][thread] = currentStepDefault;
        }
        result[i] = resultTmp[i][0];
        for (int thread = 1; thread < numThreads; thread++)
        {
            if (resultTmp[i][thread] < result[i])
                result[i] = resultTmp[i][thread];
        }
    }

    return;
}

template <class PointType>
template <class deviceType>
std::vector<float> Solver<PointType>::estimateError(const std::shared_ptr<deviceType>& deviceStatus,
                                                    float&                             maxEr)
{

    std::vector<std::vector<float>> errors(deviceStatus->GetNumberParticlesFlows());
    /*	deviceType deviceStatus2 = deviceStatus;
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
            EmptyPlaces[i].clear();
            EmptyPlaces1[i].clear();

            deviceStatus->GetFlow(i)->GetDynamicsData()->clear();
            deviceStatus2.GetFlow(i)->GetDynamicsData()->clear();

            std::vector < unsigned int > v;
            deviceStatus->GetFlow(i)->GenerateParticles(v, 1, 1, 0, deviceStatus->GetGridData(), 0,
    0);
            deviceStatus2.GetFlow(i)->GenerateParticles(v, 1, 1, 0, deviceStatus->GetGridData(), 0,
    0);

            particleGridInterface->initBoundaries(deviceStatus->GetFlow(i)->GetboundaryConditions());

            int i1 = 0;
            int i2 = deviceStatus->GetFlow(i)->GetDynamicsData()->NParticles();


            std::vector<std::vector<PointType>> v1;
            std::vector<std::vector<PointType>> v2;
            while (deviceStatus->GetFlow(i)->GetDynamicsData()->NParticles() != 0)
            {
                    EmptyPlaces1[i] =
    particleGridInterface->InCell(deviceStatus->GetFlow(i)->GetDynamicsData(),
    EmptyPlaces[i], i1, i2);

                    for (int kk = 0; kk < EmptyPlaces1[i].size(); kk++)
                    {
                            deviceStatus->GetFlow(i)->GetDynamicsData()->cellsNumbers[EmptyPlaces1[i][kk]]
    =
    deviceStatus->GetFlow(i)->GetDynamicsDataTmp()->cellsNumbers[EmptyPlaces1[i][kk]];
                    };

                    particleGridInterface->Wcalculate(deviceStatus->GetFlow(i)->GetDynamicsData(),
    W[i], i1, i2);


                    if (deviceStatus->GetFlow(i)->GetDynamicsData()->NParticles() <= 0)
                            break;

                    particleGridInterface->Grid2Particles(deviceStatus->GetFlow(i)->GetDynamicsData(),
    deviceStatus->GetGridData(), W[i], i1, i2);


                    particleMover->updateMomentums(deviceStatus->GetFlow(i)->GetDynamicsData(),
    deviceStatus->GetFlow(i)->GetAlpha(), i, i1, i2);


                    deviceStatus->GetFlow(i)->CopyDynamicsDataToTmp(i1, i2);

                    particleMover->updatePositions(deviceStatus->GetFlow(i)->GetDynamicsData(), i,
    i1, i2);

                    particleGridInterface->CheckParticlesBoundaries(EmptyPlaces[i] ,
    boundary->Getboundaries(),
    deviceStatus->GetElectrodes(), deviceStatus->GetFlow(i)->GetDynamicsDataTmp(),
    deviceStatus->GetFlow(i)->GetDynamicsData(), particleMover->GetTimeStep(i) /
    LIGHT_VELOCITY(),
    deviceStatus->GetFlow(i)->GetCharge(), deviceStatus->GetFlow(i)->GetMass(), i1, i2);

                    deviceStatus->GetFlow(i)->GetDynamicsData()->removeParticle(EmptyPlaces[i]);
                    deviceStatus2.GetFlow(i)->GetDynamicsData()->removeParticle(EmptyPlaces[i]);


                    particleMover->HalveTimeSteps();

                    int i1 = 0;
                    int i2 = deviceStatus2.GetFlow(i)->GetDynamicsData()->NParticles();

                    for (int s = 0; s < 2; s++)
                    {
                            particleGridInterface->InCell(deviceStatus2.GetFlow(i)->GetDynamicsData(),
    EmptyPlaces[i],
    i1, i2); particleGridInterface->Wcalculate(deviceStatus2.GetFlow(i)->GetDynamicsData(), W[i],
    i1, i2);
                            particleGridInterface->Grid2Particles(deviceStatus2.GetFlow(i)->GetDynamicsData(),
    deviceStatus2.GetGridData(), W[i], i1, i2);
                            particleMover->updateMomentums(deviceStatus2.GetFlow(i)->GetDynamicsData(),
    deviceStatus2.GetFlow(i)->GetAlpha(), i, i1, i2);
    deviceStatus2.GetFlow(i)->CopyDynamicsDataToTmp(i1, i2);
                            particleMover->updatePositions(deviceStatus2.GetFlow(i)->GetDynamicsData(),
    i, i1, i2);
                            particleGridInterface->CheckParticlesBoundaries(deviceStatus2.boundaries,
    deviceStatus2.conductorList, deviceStatus2.GetFlow(i)->GetDynamicsDataTmp(),
    deviceStatus2.GetFlow(i)->GetDynamicsData(), particleMover->GetTimeStep(i) /
    LIGHT_VELOCITY(),
    deviceStatus2.GetFlow(i)->GetCharge(), deviceStatus2.GetFlow(i)->GetMass(), i1, i2);
                    }
                    particleMover->DoubleTimeSteps();


                    v1 = deviceStatus->GetFlow(i)->GetDynamicsData()->GetDataVector();
                    v2 = deviceStatus2.GetFlow(i)->GetDynamicsData()->GetDataVector();


                    PointType maxTol = -1;
                    for (int i = 0; i < v1.size(); i++)
                    {
                            for (int j = 0; j < v1[i].size(); j++)
                            {
                                    PointType tol = std::abs(v1[i][j] - v2[i][j]) / 3;
                                    if (tol > maxTol)
                                            maxTol = tol;
                            }

                    };
                    if (maxTol>0)
                            errors[i].push_back(maxTol);
            }



    }


    std::vector<float> result;
    maxEr=-1;
    float maxErLoc = -1;

    /*for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
            maxErLoc = -1;
            for (int j = 0; j < errors[i].size() - 1; j++)
            {
                    if (errors[i][j]>maxErLoc)
                            maxErLoc = errors[i][j];
            }
            if (maxEr > maxErLoc)
                    maxEr = maxErLoc;
            result.push_back(maxErLoc);
    }*/

    /*	for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
            {
                    for (int j = 0; j < errors[i].size() - 1; j++)
                    {
                            if (errors[i][j]>maxEr)
                                    maxEr = errors[i][j];
                    }
            }*/
    return errors[0];
}

template <class PointType>
template <class deviceType>
void Solver<PointType>::ErrorEstimateEvent(const std::shared_ptr<deviceType>& deviceStatus)
{

    PointType Lmax = 0;
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        PointType l = deviceStatus->GetFlow(i)->GetEmitterDevice()->GetSourceSize();
        if (l > Lmax)
            Lmax = l;
    }
    double progressLoc;

    InitFieldSolver(deviceStatus, progressLoc);
    int memorySize = 1e5;

    std::vector<int> boundaryPoints = fieldSolver->getNearBoundarypoints_all();

    particleGridInterface->init(deviceStatus->GetNumberParticlesFlows(),
                                deviceStatus->GetGridData(), deviceStatus->Getmesh()->templNumb,
                                deviceStatus->Getmesh()->flagMatrix, boundaryPoints,
                                deviceStatus->GetDomainBoundary(), memorySize, 1, 10000);

    fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                               deviceStatus->GetboundaryConditions(), 0, progressLoc);

    W.resize(deviceStatus->GetNumberParticlesFlows());
    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
        W[i]   = new PointType[memorySize][9];

    Wtmp.resize(deviceStatus->GetNumberParticlesFlows());
    for (int i  = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
        Wtmp[i] = new PointType[memorySize][9];

    EmptyPlaces.resize(deviceStatus->GetNumberParticlesFlows());
    EmptyPlaces1.resize(deviceStatus->GetNumberParticlesFlows());

    std::vector<int> flag;

    // emissionCurrentSolver->reset();

    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
    {
        deviceStatus->GetElectrodes()[i]->ResetPower();
    }
    int iteration = 0;

    double dH = deviceStatus->GetGridData()->GetMaxSixe();

    std::vector<PointType> alpas(deviceStatus->GetNumberParticlesFlows());

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
    {
        flag.push_back(1);
        deviceStatus->GetFlow(i)->ReserveMemory(memorySize);
        std::vector<unsigned int> t_empty;
        deviceStatus->GetFlow(i)->GenerateParticles(t_empty, 1, 1, 0, deviceStatus->GetGridData(),
                                                    0, 0);
        particleGridInterface->SearchStartCells(i, deviceStatus->GetFlow(i)->GetDynamicsData());
        std::vector<unsigned int> realCells =
            deviceStatus->GetFlow(i)->GetDynamicsData()->Get_startCellNumbers();
        //	emissionCurrentSolver->init(deviceStatus->GetGridData(), dH,
        // deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(),
        // deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSource(),
        // &particleGridInterface,
        // deviceStatus->GetGridData()->GetType(), realCells, Lmax);
        deviceStatus->GetFlow(i)->GetDynamicsData()->clear();
        std::vector<unsigned int> v;
        alpas[i] = deviceStatus->GetFlow(i)->GetAlpha();
    }

    particleMover->init(alpas, 1);

    float maxEr;
    errors = estimateError(deviceStatus, maxEr);
}


template class Solver<float>;
template class Solver<double>;

template void Solver<double>::TimeStepEstimate<device2daxsdouble>(
        std::vector<double>& result, const std::shared_ptr<device2daxsdouble>& deviceStatus,
        int numThreads, double& progress);

template void Solver<float>::TimeStepEstimate<device2daxsfloat>(
        std::vector<double>& result, const std::shared_ptr<device2daxsfloat>& deviceStatus,
        int numThreads, double& progress);

template void Solver<double>::TimeStepEstimate<device2ddouble>(
        std::vector<double>& result, const std::shared_ptr<device2ddouble>& deviceStatus,
        int numThreads, double& progress);

template void
Solver<float>::TimeStepEstimate<device2dfloat>(std::vector<double>&                  result,
                                               const std::shared_ptr<device2dfloat>& deviceStatus,
                                               int numThreads, double& progress);

template void Solver<double>::TimeStepEstimate<device2dpolardouble>(
        std::vector<double>& result, const std::shared_ptr<device2dpolardouble>& deviceStatus,
        int numThreads, double& progress);

template void Solver<float>::TimeStepEstimate<device2dpolarfloat>(
        std::vector<double>& result, const std::shared_ptr<device2dpolarfloat>& deviceStatus,
        int numThreads, double& progress);

template void Solver<double>::TimeStepEstimate<device3dExtrdouble>(
        std::vector<double>& result, const std::shared_ptr<device3dExtrdouble>& deviceStatus,
        int numThreads, double& progress);

template void Solver<float>::TimeStepEstimate<device3dExtrfloat>(
        std::vector<double>& result, const std::shared_ptr<device3dExtrfloat>& deviceStatus,
        int numThreads, double& progress);

template void Solver<double>::TimeStepEstimate<device3ddouble>(
        std::vector<double>& result, const std::shared_ptr<device3ddouble>& deviceStatus,
        int numThreads, double& progress);

template void Solver<double>::ErrorEstimateEvent<device2daxsdouble>(
        const std::shared_ptr<device2daxsdouble>& deviceStatus);

template void Solver<float>::ErrorEstimateEvent<device2daxsfloat>(
        const std::shared_ptr<device2daxsfloat>& deviceStatus);

template void Solver<double>::ErrorEstimateEvent<device2ddouble>(
        const std::shared_ptr<device2ddouble>& deviceStatus);

template void Solver<float>::ErrorEstimateEvent<device2dfloat>(
        const std::shared_ptr<device2dfloat>& deviceStatus);

template void Solver<double>::ErrorEstimateEvent<device2dpolardouble>(
        const std::shared_ptr<device2dpolardouble>& deviceStatus);

template void Solver<float>::ErrorEstimateEvent<device2dpolarfloat>(
        const std::shared_ptr<device2dpolarfloat>& deviceStatus);

template void Solver<double>::ErrorEstimateEvent<device3dExtrdouble>(
        const std::shared_ptr<device3dExtrdouble>& deviceStatus);

template void Solver<float>::ErrorEstimateEvent<device3dExtrfloat>(
        const std::shared_ptr<device3dExtrfloat>& deviceStatus);

template void Solver<double>::ErrorEstimateEvent<device3ddouble>(
        const std::shared_ptr<device3ddouble>& deviceStatus);

template std::vector<float> Solver<double>::estimateError<device2daxsdouble>(
        const std::shared_ptr<device2daxsdouble>& deviceStatus, float& maxEr);

template std::vector<float> Solver<float>::estimateError<device2daxsfloat>(
        const std::shared_ptr<device2daxsfloat>& deviceStatus, float& maxEr);

template std::vector<float>
Solver<double>::estimateError<device2ddouble>(const std::shared_ptr<device2ddouble>& deviceStatus,
                                              float&                                 maxEr);

template std::vector<float>
Solver<float>::estimateError<device2dfloat>(const std::shared_ptr<device2dfloat>& deviceStatus,
                                            float&                                maxEr);

template std::vector<float> Solver<double>::estimateError<device2dpolardouble>(
        const std::shared_ptr<device2dpolardouble>& deviceStatus, float& maxEr);

template std::vector<float> Solver<float>::estimateError<device2dpolarfloat>(
        const std::shared_ptr<device2dpolarfloat>& deviceStatus, float& maxEr);

template std::vector<float> Solver<double>::estimateError<device3dExtrdouble>(
        const std::shared_ptr<device3dExtrdouble>& deviceStatus, float& maxEr);

template std::vector<float> Solver<float>::estimateError<device3dExtrfloat>(
        const std::shared_ptr<device3dExtrfloat>& deviceStatus, float& maxEr);

template std::vector<float>
Solver<double>::estimateError<device3ddouble>(const std::shared_ptr<device3ddouble>& deviceStatus,
                                              float&                                 maxEr);