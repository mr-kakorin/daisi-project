#include "ParticlesFlow.h"
#include "BoundaryConditions.h"
#include "Dmath.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice3d.h"
#include "Geom.h"
#include "Particle.h"
#include "ParticleSource.h"
#include "Tools.h"
#include <cmath>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <common_tools/constants.h>

template class ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>;
template class ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>;

template class ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>;
template class ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>;

template class ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>;
template class ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>;

template class ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>;
template class ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>;

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetEmittanceData(
    std::vector<std::vector<float>>& data, int flagNumber, int emFlag)
{
    Emittances[0][flagNumber]->GetEmittanceData(data, emFlag, mass,
                                                EmitterDeviceStatus->GetLambda());
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::AddEmittancesData(int thread)
{
    for (int i = 0; i < Emittances[thread].size(); i++)
        Emittances[thread][i]->InsertParticelsEmittances(
            DynamicsDataParallel[thread].get(), writeIndexes[thread][i], EmittancesList[i]);

    for (int i = 0; i < Emittances[thread].size(); i++)
        writeIndexes[thread][i].clear();
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
bool ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::isConfigure()
{
    if (DistributionStyle < 5 &&
        !EmitterDeviceStatus->GetParticleSources()[0]->sourceSurface.size())
        return false;
    return true;
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::MergeEmittancesData()
{
    for (int thread = 1; thread < Emittances.size(); thread++)
    {
        for (int i = 0; i < Emittances[0].size(); i++)
        {
            Emittances[0][i]->InsertParticels(Emittances[thread][i].get());
            Emittances[thread][i]->clear();
        }
    }
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetRmsEmittances(
    std::vector<float>& emittance)
{
    std::vector<std::vector<PointType>> meas(4);
    std::vector<PointType>              pos;

    for (int i = 0; i < 4; i++)
        meas[i].resize(5);

    emittance.resize(3);

    float  n      = 0;
    float  zAv    = 0;
    double betaAv = 0;

    for (int i = 0; i < DynamicsDataParallel.size(); i++)
    {
        DynamicsDataParallel[i]->GetAveragePositions(pos);
        zAv    = zAv + pos[0];
        betaAv = betaAv + pos[1];
        n      = n + DynamicsDataParallel[i]->NParticles();
    }
    betaAv = betaAv / DynamicsDataParallel.size();
    zAv    = zAv / DynamicsDataParallel.size();

    for (int i = 0; i < DynamicsDataParallel.size(); i++)
    {
        DynamicsDataParallel[i]->GetBeamMeasuriments(meas, mass, zAv, betaAv);
    };

    for (int i = 0; i < 3; i++)
    {
        float xx     = meas[i][0] / n - (meas[i][2] / n) * (meas[i][2] / n);
        float pp     = meas[i][1] / n - (meas[i][3] / n) * (meas[i][3] / n);
        float xxpp   = meas[i][4] / n - meas[i][2] * meas[i][3] / (n * n);
        emittance[i] = sqrt(xx * pp - xxpp * xxpp);
    };
    emittance[0] = 4 * betaAv * emittance[0] * 1e5;
    emittance[1] = 4 * betaAv * emittance[1] * 1e5;
    emittance[2] = 4 * 1e9 * emittance[2] / (betaAv * commtools::LIGHT_VELOCITY());

    for (int i = 0; i < 3; i++)
    {
        if (std::isinf(emittance[i]) || std::isnan(emittance[i]) || emittance[i] < 0)
            emittance[i] = 0;
    }
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<double>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetEmittancesList()
{
    return EmittancesList;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::AddEmittance(double param)
{
    EmittancesList.push_back(param);
    Emittances[0].push_back(std::shared_ptr<ParticlesDataType>(new ParticlesDataType()));
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::checkEmittances(int thread,
                                                                                     int i1, int i2)
{

    for (int i = 0; i < Emittances[thread].size(); i++)
    {
        std::vector<int> numbers = DynamicsDataParallel[thread]->CheckEmittanceCondition(
            EmittancesList[i], DynamicsDataParallelTmp[thread].get(), i1, i2);
        writeIndexes[thread][i].insert(writeIndexes[thread][i].end(), numbers.begin(),
                                       numbers.end());
    }
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
long long ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetMemorySize()
{
    volatile long long result = 0;
    for (int i = 0; i < DynamicsDataParallel.size(); i++)
        result = result + DynamicsDataParallel[i]->GetMemorySize();

    for (int i = 0; i < DynamicsDataParallelTmp.size(); i++)
        result = result + DynamicsDataParallelTmp[i]->GetMemorySize();

    result = result + DynamicsData->GetMemorySize();
    result = result + tmpData->GetMemorySize();
    result = result + newParticlesZ0->GetMemorySize();

    return result;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
int ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetNumberOfThreads()
{
    return DynamicsDataParallel.size();
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
int ParticlesFlow<ParticlesDataType, EmitterDeviceType,
                  PointType>::GetNumberOfParticles() // total number of particles in
                                                     // DynamicsDataParallel
{
    int res = 0;
    for (int i = 0; i < DynamicsDataParallel.size(); i++)
        res = res + DynamicsDataParallel[i]->NParticles();
    return res;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GenerateSyncParticle()
{
    EmitterDeviceStatus->GenerateSyncParticle(DynamicsData, mass);
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::CalculateFlowCurrent()
{
    double res = 0;

    for (int i = 0; i < DynamicsDataParallel.size(); i++)
        res = res + std::abs(DynamicsDataParallel[i]->GetTotalCurrent());

    EmitterDeviceStatus->SetFlowCurrent(res);
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
int ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetMaxParticles()
{
    return EmitterDeviceStatus->GetMaxParticles();
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::DivideForThreads(
    int numThreads){

};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::clearParallel()
{
    for (int j = 0; j < DynamicsDataParallel.size(); j++)
        DynamicsDataParallel[j]->clear();
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::InitParallel(int size,
                                                                                  int numThreads,
                                                                                  int flagMethod)
{
    newParticlesZ0 = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));

    Emittances.resize(numThreads);

    for (int j = 0; j < numThreads; j++)
    {
        Emittances[j].resize(Emittances[0].size());

        for (int i = 0; i < Emittances[0].size(); i++)
        {
            Emittances[j][i] =
                std::shared_ptr<ParticlesDataType>(new ParticlesDataType(flagMethod));
            Emittances[j][i]->clear();
        }
    };

    int nowThreads = DynamicsDataParallel.size();

    if (nowThreads < numThreads)
    {
        DynamicsDataParallel.resize(numThreads);
        for (int j = nowThreads; j < numThreads; j++)
        {
            DynamicsDataParallel[j] =
                std::shared_ptr<ParticlesDataType>(new ParticlesDataType(flagMethod));
            DynamicsDataParallel[j]->ReserveMemory(size / numThreads);
        }
        int nParticlesPerThread = GetNumberOfParticles() / numThreads;
        for (int j = 0; j < nowThreads; j++)
        {
            int rebalance = (DynamicsDataParallel[j]->NParticles() - nParticlesPerThread) /
                            (numThreads - nowThreads);
            for (int i = nowThreads; i < numThreads - 1; i++)
                DynamicsDataParallel[i]->RecombinateParticles(DynamicsDataParallel[j].get(),
                                                              rebalance);

            DynamicsDataParallel[numThreads - 1]->RecombinateParticles(
                DynamicsDataParallel[j].get(), DynamicsDataParallel[j]->NParticles());
        }
    };

    if (nowThreads > numThreads)
    {
        for (int j = numThreads; j < nowThreads; j++)
        {
            int rebalance = (DynamicsDataParallel[j]->NParticles()) / (numThreads);
            for (int i = 0; i < numThreads - 1; i++)
                DynamicsDataParallel[i]->RecombinateParticles(DynamicsDataParallel[j].get(),
                                                              rebalance);

            DynamicsDataParallel[numThreads - 1]->RecombinateParticles(
                DynamicsDataParallel[j].get(), DynamicsDataParallel[j]->NParticles());
        };
        DynamicsDataParallel.resize(numThreads);
    };

    writeIndexes.resize(numThreads);

    indesexPerThread.clear();

    indesexPerThread.resize(numThreads);

    int n;
    if (DistributionStyle < 5)
        n = EmitterDeviceStatus->GetnParticlesEmitter();
    else
        n = EmitterDeviceStatus->GetnParticlesBunch();

    int perThread     = n / numThreads;
    int perLastThread = n - perThread * (numThreads - 1);

    for (int i = 0; i < numThreads; i++)
    {
        writeIndexes[i].resize(Emittances[i].size());
    }

    std::vector<int> allIndexes;

    for (int i = 0; i < n; i++)
        allIndexes.push_back(i);

    for (int i = 0; i < numThreads; i++)
    {
        int nLoc = perThread;
        if (i == numThreads - 1)
            nLoc = perLastThread;
        for (int j = 0; j < nLoc; j++)
        {
            int index = rand() % allIndexes.size();
            indesexPerThread[i].push_back(allIndexes[index]);
            allIndexes.erase(allIndexes.begin() + index);
        };
    }
    for (int i = 0; i < numThreads; i++)
    {
        std::sort(indesexPerThread[i].begin(), indesexPerThread[i].end());
    };
}

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::InitParallelTmpBlocks(
    int numThreads, int blockSize, int solverType)
{
    DynamicsDataParallelTmp.resize(numThreads);

    for (int i = 0; i < numThreads; i++)
    {
        DynamicsDataParallelTmp[i] =
            std::shared_ptr<ParticlesDataType>(new ParticlesDataType(solverType));
        DynamicsDataParallelTmp[i]->ReserveMemory(blockSize * 2);
        DynamicsDataParallelTmp[i]->resize(blockSize * 2);
    }
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<double>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetFlowMCNumbers()
{
    std::vector<double> result(3);
    result[0] = mass / commtools::PROTON_MASS();
    result[1] = charge / std::abs(commtools::ELECTRON_CHARGE());
    result[2] = maxTime;
    return result;
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::SetFlowMCNumbers(
    std::vector<double> numbers)
{
    mass    = commtools::PROTON_MASS() * numbers[0];
    charge  = std::abs(commtools::ELECTRON_CHARGE()) * numbers[1];
    maxTime = numbers[2];
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
int ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetParticlesType()
{
    return ParticlesType;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
int ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDistributionStyle()
{
    return DistributionStyle;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
PointType ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetMass()
{
    return mass;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
PointType ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetCharge()
{
    return charge;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
PointType ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetAlpha()
{
    return charge / (mass * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY());
};
// EmitterDeviceType*  EmitterDeviceStatus;

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<EmitterDeviceType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetEmitterDevice()
{
    return EmitterDeviceStatus;
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<BoundaryConditions>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetboundaryConditions()
{
    return boundaryConditions;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<ParticlesDataType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsData()
{
    return DynamicsData;
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<ParticlesDataType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsDataStart()
{
    if (DistributionStyle == 5)
        return newParticlesZ0;
    else
        return DynamicsData;
}

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<ParticlesDataType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsDataTmp()
{
    return tmpData;
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<ParticlesDataType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsData(int thread)
{
    return DynamicsDataParallel[thread];
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::shared_ptr<ParticlesDataType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsDataTmp(int thread)
{
    return DynamicsDataParallelTmp[thread];
}

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<std::shared_ptr<ParticlesDataType>>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetDynamicsDataParallelArray()
{
    return DynamicsDataParallel;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<float>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetData(int Index)
{
    return DynamicsData->GetData(Index, mass, charge);
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<void*> ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetData()
{
    return DynamicsData->GetData();
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::ParticlesFlow()
{
    EmitterDeviceStatus = std::shared_ptr<EmitterDeviceType>(new EmitterDeviceType());
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::ReserveMemory(int size){
    //	DynamicsData->ReserveMemory()(size);
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::ParticlesFlow(
    int ParticlesTypeIn, int DistributionStyleIn, double MassNumber, double ChargeNumber,
    std::vector<int> boundariesList) //����������� ��� �����
{
    Emittances.resize(1);
    EmitterDeviceStatus =
        std::shared_ptr<EmitterDeviceType>(new EmitterDeviceType(DistributionStyleIn));
    ParticlesType      = ParticlesTypeIn;
    DistributionStyle  = DistributionStyleIn;
    mass               = commtools::PROTON_MASS() * MassNumber;
    charge             = std::abs(commtools::ELECTRON_CHARGE()) * ChargeNumber;
    maxTime            = 200;
    boundaryConditions = std::shared_ptr<BoundaryConditions>(new BoundaryConditions());
    boundaryConditions->SetDefaultConditionsList(boundariesList);

    newParticlesZ0 = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
    tmpData        = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
    DynamicsData   = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
bool ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::CheckTimeLimit(int thread)
{
    if (maxTime * 1e-9 * commtools::LIGHT_VELOCITY() < DynamicsDataParallel[thread]->Time)
        return true;
    return false;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::ParticlesFlow(
    int ParticlesTypeIn, int DistributionStyleIn,
    std::vector<int> boundariesList) //����������� ����������
{
    Emittances.resize(1);
    EmitterDeviceStatus = std::shared_ptr<EmitterDeviceType>(new EmitterDeviceType());
    ParticlesType       = ParticlesTypeIn;
    DistributionStyle   = DistributionStyleIn;
    mass                = commtools::ELECTRON_MASS();
    charge              = commtools::ELECTRON_CHARGE();
    maxTime             = 200;
    boundaryConditions  = std::shared_ptr<BoundaryConditions>(new BoundaryConditions());
    boundaryConditions->SetDefaultConditionsList(boundariesList);

    newParticlesZ0 = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
    tmpData        = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
    DynamicsData   = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<double>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetFlowProperties()
{
    std::vector<double> properties;
    properties.push_back(double(ParticlesType));
    properties.push_back(double(DistributionStyle));
    properties.push_back(mass);
    properties.push_back(charge);
    properties.push_back(maxTime);
    return properties;
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
std::vector<unsigned int>
ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GetNewIndexes(
    std::vector<unsigned int> EmptyPlaces)
{
    int particlesNumber = EmitterDeviceStatus->GetNumbersOfParticlesGeneration();
    std::vector<unsigned int> result;

    for (int i = 0; i < DynamicsData->NParticles(); i++)
    {
        if (int(DynamicsData->flagEmitted[i]) != 1)
            result.push_back(i);
    }
    /*if (particlesNumber<EmptyPlaces.size())
    {
    for (int i = 0; i < particlesNumber; i++)
    result[i] = EmptyPlaces[i];
    return result;
    }

    int k = 0;

    for (int i = 0; i < EmptyPlaces.size(); i++)
    result[i] = EmptyPlaces[i];

    for (int i = int(EmptyPlaces.size()); i < particlesNumber; i++)
    {
    result[i] = DynamicsData->NParticles() - (particlesNumber - EmptyPlaces.size()) + k;
    k++;
    }*/
    return result;
};

template void
ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>::GenerateParticles<
    GridData2daxs<float>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                          int stepNumber, const std::shared_ptr<GridData2daxs<float>>& grid,
                          int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>::GenerateParticles<
    GridData2daxs<double>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridData2daxs<double>>& grid,
                           int flagLocate, int flagDistrib);

template void ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>::GenerateParticles<
    GridData2d<float>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                       int stepNumber, const std::shared_ptr<GridData2d<float>>& grid,
                       int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>::GenerateParticles<
    GridData2d<double>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                        int stepNumber, const std::shared_ptr<GridData2d<double>>& grid,
                        int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>::GenerateParticles<
    GridData2dpolar<float>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                            int stepNumber, const std::shared_ptr<GridData2dpolar<float>>& grid,
                            int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>::GenerateParticles<
    GridData2dpolar<double>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                             int stepNumber, const std::shared_ptr<GridData2dpolar<double>>& grid,
                             int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>::GenerateParticles<
    GridData3d<double>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                        int stepNumber, const std::shared_ptr<GridData3d<double>>& grid,
                        int flagLocate, int flagDistrib);

template void ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>::GenerateParticles<
    GridData3d<float>>(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                       int stepNumber, const std::shared_ptr<GridData3d<float>>& grid,
                       int flagLocate, int flagDistrib);


template <class ParticlesDataType, class EmitterDeviceType, class PointType>
template <class GridDataType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GenerateParticles(
    std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt, int stepNumber,
    const std::shared_ptr<GridDataType>& grid, int flagLocate, int flagDistrib)
{

    if (DistributionStyle == 5)
        EmitterDeviceStatus->PreliminaryGeneration(newParticlesZ0,
                                                   mass); // ellipse generation and soft
    else
    {
        std::vector<int> allIndexes;
        int              n = EmitterDeviceStatus->GetnParticlesEmitter();
        for (int i = 0; i < n; i++)
            allIndexes.push_back(i);
        EmitterDeviceStatus->GenerateParticles(allIndexes, 0, 1, EmptyPlaces, DynamicsData, mass,
                                               charge, flagClear, dt, stepNumber, grid, flagLocate,
                                               flagDistrib);
    };
};

template void
ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>::GenerateParticlesThreaded<
    GridData2daxs<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                          int flagClear, double dt, int stepNumber,
                          const std::shared_ptr<GridData2daxs<float>>& grid, int flagLocate,
                          int flagDistrib);

template void ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>::
    GenerateParticlesThreaded<GridData2daxs<double>>(
        int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear,
        double dt, int stepNumber, const std::shared_ptr<GridData2daxs<double>>& grid,
        int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>::GenerateParticlesThreaded<
    GridData2d<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                       int flagClear, double dt, int stepNumber,
                       const std::shared_ptr<GridData2d<float>>& grid, int flagLocate,
                       int flagDistrib);

template void
ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>::GenerateParticlesThreaded<
    GridData2d<double>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                        int flagClear, double dt, int stepNumber,
                        const std::shared_ptr<GridData2d<double>>& grid, int flagLocate,
                        int flagDistrib);

template void
ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>::GenerateParticlesThreaded<
    GridData2dpolar<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                            int flagClear, double dt, int stepNumber,
                            const std::shared_ptr<GridData2dpolar<float>>& grid, int flagLocate,
                            int flagDistrib);

template void
ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>::GenerateParticlesThreaded<
    GridData2dpolar<double>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                             int flagClear, double dt, int stepNumber,
                             const std::shared_ptr<GridData2dpolar<double>>& grid, int flagLocate,
                             int flagDistrib);

template void
ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>::GenerateParticlesThreaded<
    GridData3d<double>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                        int flagClear, double dt, int stepNumber,
                        const std::shared_ptr<GridData3d<double>>& grid, int flagLocate,
                        int flagDistrib);

template void
ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>::GenerateParticlesThreaded<
    GridData3d<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                       int flagClear, double dt, int stepNumber,
                       const std::shared_ptr<GridData3d<float>>& grid, int flagLocate,
                       int flagDistrib);

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
template <class GridDataType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GenerateParticlesThreaded(
    int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
    int stepNumber, const std::shared_ptr<GridDataType>& grid, int flagLocate, int flagDistrib)
{
    if (DistributionStyle == 5)
        EmitterDeviceStatus->GenerateParticlesLinac(
            0, thread, numThreads, EmptyPlaces, DynamicsDataParallel[thread], mass,
            Dmath::sign(charge), flagClear, dt, stepNumber, grid, flagLocate, flagDistrib);
    else
        EmitterDeviceStatus->GenerateParticles(
            indesexPerThread[thread], thread, numThreads, EmptyPlaces, DynamicsDataParallel[thread],
            mass, charge, flagClear, dt, stepNumber, grid, flagLocate, flagDistrib);
};

template void ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>::
    GenerateParticlesThreadedTest<GridData2daxs<float>>(
        int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear,
        double dt, int stepNumber, const std::shared_ptr<GridData2daxs<float>>& grid,
        int flagLocate, int flagDistrib);

template void ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>::
    GenerateParticlesThreadedTest<GridData2daxs<double>>(
        int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear,
        double dt, int stepNumber, const std::shared_ptr<GridData2daxs<double>>& grid,
        int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>::GenerateParticlesThreadedTest<
    GridData2d<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                       int flagClear, double dt, int stepNumber,
                       const std::shared_ptr<GridData2d<float>>& grid, int flagLocate,
                       int flagDistrib);

template void
ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>::GenerateParticlesThreadedTest<
    GridData2d<double>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                        int flagClear, double dt, int stepNumber,
                        const std::shared_ptr<GridData2d<double>>& grid, int flagLocate,
                        int flagDistrib);

template void ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>::
    GenerateParticlesThreadedTest<GridData2dpolar<float>>(
        int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear,
        double dt, int stepNumber, const std::shared_ptr<GridData2dpolar<float>>& grid,
        int flagLocate, int flagDistrib);

template void ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>::
    GenerateParticlesThreadedTest<GridData2dpolar<double>>(
        int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear,
        double dt, int stepNumber, const std::shared_ptr<GridData2dpolar<double>>& grid,
        int flagLocate, int flagDistrib);

template void
ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>::GenerateParticlesThreadedTest<
    GridData3d<double>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                        int flagClear, double dt, int stepNumber,
                        const std::shared_ptr<GridData3d<double>>& grid, int flagLocate,
                        int flagDistrib);

template void
ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>::GenerateParticlesThreadedTest<
    GridData3d<float>>(int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces,
                       int flagClear, double dt, int stepNumber,
                       const std::shared_ptr<GridData3d<float>>& grid, int flagLocate,
                       int flagDistrib);

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
template <class GridDataType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::GenerateParticlesThreadedTest(
    int thread, int numThreads, std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
    int stepNumber, const std::shared_ptr<GridDataType>& grid, int flagLocate, int flagDistrib)
{
    if (DistributionStyle == 5)
        EmitterDeviceStatus->GenerateParticlesLinac(
            1, thread, numThreads, EmptyPlaces, DynamicsDataParallel[thread], mass,
            Dmath::sign(charge), flagClear, dt, stepNumber, grid, flagLocate, flagDistrib);
    else if (stepNumber == 0)
        EmitterDeviceStatus->GenerateParticles(
            indesexPerThread[thread], thread, numThreads, EmptyPlaces, DynamicsDataParallel[thread],
            mass, charge, flagClear, dt, stepNumber, grid, flagLocate, flagDistrib);
};

/*template <class GridDataType>
void GenerateParticlesVector(double dt, int stepNumber)
{
newParticles = ParticlesDataType();
EmitterDeviceStatus->GenerateParticles(std::vector <unsigned int > {}, newParticles, mass,
Dmath::sign(charge), 1, dt,
stepNumber);
};*/
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::setNewParticles(
    const std::vector<unsigned int>& EmptyPlaces)
{
    // DynamicsData->setNewParticles(EmptyPlaces, newParticles.get());
}
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::CopyDynamicsDataToTmp(int i1,
                                                                                           int i2)
{
    DynamicsData->FastCopy(tmpData.get(), i1, i2);
};
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::CopyDynamicsDataToTmpThreaded(
    int thread, int i1, int i2)
{
    DynamicsDataParallel[thread]->FastCopy(DynamicsDataParallelTmp[thread].get(), i1, i2);
};

template void ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>::serialize<
    boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                     const unsigned int               file_version);

template void ParticlesFlow<Particles3dcil<float>, EmitterDevice2daxs<float>, float>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3dcil<double>, EmitterDevice2daxs<double>, double>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2d<float>, EmitterDevice2d<float>, float>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2d<double>, EmitterDevice2d<double>, double>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2dpolar<float>, EmitterDevice2d<float>, float>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles2dpolar<double>, EmitterDevice2d<double>, double>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3d<double>, EmitterDevice3d<double>, double>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);
template void ParticlesFlow<Particles3d<float>, EmitterDevice3d<float>, float>::serialize<
    boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                     const unsigned int               file_version);

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
template <class Archive>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::save(Archive& ar,
                                                                          const unsigned int) const
{
    ar& ParticlesType;
    ar& mass;
    ar& charge;
    ar& EmitterDeviceStatus;
    //	ar & DynamicsData;
    ar& boundaryConditions;
    ar& DistributionStyle;
    ar& newParticlesZ0;
    ar& maxTime;
    ar& EmittancesList;
    ar& Emittances[0];
};

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
template <class Archive>
void ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>::load(Archive& ar,
                                                                          const unsigned int)
{
    ar& ParticlesType;
    ar& mass;
    ar& charge;
    ar& EmitterDeviceStatus;
    //	ar & DynamicsData;
    ar& boundaryConditions;
    ar& DistributionStyle;
    ar& newParticlesZ0;
    ar& maxTime;
    ar& EmittancesList;
    Emittances.resize(1);
    ar& Emittances[0];
    DynamicsData = std::shared_ptr<ParticlesDataType>(new ParticlesDataType(0));
};