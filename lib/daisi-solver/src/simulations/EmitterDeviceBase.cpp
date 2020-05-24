#include "EmitterDeviceBase.h"
#include "BoundaryContainer2d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "Geom.h"
#include "ParticleSource.h"
#include "Tools.h"
#include <Constants.h>

static std::default_random_engine             generator;
static std::normal_distribution<double>       distribution(0, 1);
static std::uniform_real_distribution<double> distributionUn(0, 1);

template class EmitterDeviceBase<float>;
template class EmitterDeviceBase<double>;

template <class PointType>
void EmitterDeviceBase<PointType>::GetSliceIndexes(std::vector<int>& sliceIndexesParallel,
                                                   int flagTest, double t1, double t2,
                                                   double& phMin, double& phMax, int numThreads,
                                                   int thread)
{

    std::vector<int> sliceIndexes;

    double frequency = DistribParams[3];

    double z1 = 2 * PI() * t1 * frequency;
    double z2 = 2 * PI() * t2 * frequency;

    int nBunches = NumbersParams[1];

    if (z2 < nBunches * 2 * PI())
    {
        while (z1 > 2 * PI())
            z1 = z1 - 2 * PI();

        while (z2 > 2 * PI())
            z2 = z2 - 2 * PI();
    }
    // double z1 = fmod(t1*frequency, 2 * PI()) - PI();
    // double z2 = fmod(t2*frequency, 2 * PI()) - PI();
    for (int i = 0; i < Data[4].size(); i++)
    {
        if (Data[4][i] > z1 && Data[4][i] < z2)
            sliceIndexes.push_back(i);
    };
    /*
    int i = 0;
    while (Data[4][i]>z2)
    {
            i++;
            if (i >= Data[4].size() - 1)
                    break;
    }
    while (Data[4][i]>z1 && Data[4][i] < z2)
    {
            sliceIndexes.push_back(i);
            i++;
            if (i >= Data[4].size() - 1)
                    break;
    }
    */
    int s   = sliceIndexes.size() / nBunches;
    int rem = sliceIndexes.size() - s;

    for (int i = 0; i < rem; i++)
    {
        int r = rand() % sliceIndexes.size();
        sliceIndexes.erase(sliceIndexes.begin() + r);
    };

    phMin = 100000;
    phMax = -100000;

    for (int i = 0; i < sliceIndexes.size(); i++)
    {
        if (Data[4][sliceIndexes[i]] > phMax)
            phMax = Data[4][sliceIndexes[i]];

        if (Data[4][sliceIndexes[i]] < phMin)
            phMin = Data[4][sliceIndexes[i]];
    };

    // if (Data[4][particlesPerBunch]>z1)//there is no bunch coming at this time, just remove
    // Emptyplaces
    //{//emitpiriod?
    //}

    int perThread = sliceIndexes.size() / numThreads;

    int j1 = thread * perThread;
    int j2 = (thread + 1) * perThread;
    if (thread == numThreads - 1)
    {
        j2 = sliceIndexes.size();
    };

    sliceIndexesParallel.resize(j2 - j1);

    for (int j = j1; j < j2; j++)
    {
        sliceIndexesParallel[j - j1] = sliceIndexes[j];
    };

    int totalParticles;

    if (flagTest == 1)
        totalParticles = std::min(int(sliceIndexesParallel.size()), 2);
    else
        totalParticles = sliceIndexesParallel.size();

    sliceIndexesParallel.resize(totalParticles);
};

template <class PointType>
void EmitterDeviceBase<PointType>::GenerateEllipses(double restMass)
{
    int particlesPerBunch = NumbersParams[0] * NumbersParams[1];

    double energyAverage = DistribParams[1];
    double restEnergy    = -restMass * LIGHT_VELOCITY() * LIGHT_VELOCITY() /
                        ELECTRON_CHARGE(); // ������� ����� ������� � ���������������

    double gamma = (restEnergy + std::abs(energyAverage)) / restEnergy;
    double beta  = sqrt(gamma * gamma - 1) / gamma;

    double P = gamma * beta;

    // Xparam[0] = (DistribParams[5] - DistribParams[4]) / 2; //phaseSpread

    // Yparam[0] = P*DistribParams[2] / 100;//energySpread

    double emittanceXunnorm = DistribParams[6] / beta;
    double emittanceYunnorm = DistribParams[9] / beta;

    double Xparam[3];
    double Yparam[3];
    double Angles[3];
    double e[2];

    double alpha[3];

    Xparam[0] = DistribParams[7];        // double X
    Yparam[0] = DistribParams[8];        // dX
    e[0]      = emittanceXunnorm * 1e-6; // em

    Xparam[1] = DistribParams[10];       // double Y
    Yparam[1] = DistribParams[11];       // dY
    e[1]      = emittanceYunnorm * 1e-6; // em

    double gammaE[2];
    double betaE[2];
    double phi[2];
    double k[2];

    double Major[2];
    double Minor[2];

    for (int i = 0; i < 2; i++)
    {
        //	e[i] = Xparam[i] * Yparam[i] / sqrt(1 + alpha[i] * alpha[i]);

        alpha[i]  = sqrt(Xparam[i] * Xparam[i] * Yparam[i] * Yparam[i] - e[i] * e[i]) / e[i];
        gammaE[i] = Yparam[i] * Yparam[i] / std::abs(e[i]);
        betaE[i]  = Xparam[i] * Xparam[i] / std::abs(e[i]);

        phi[i] = atan(2 * alpha[i] / (gammaE[i] - betaE[i])) / 2;

        k[i] = tan(phi[i]);

        double x = sqrt(std::abs(e[i]) / (gammaE[i] + 2 * alpha[i] * k[i] + betaE[i] * k[i] * k[i]));
        double y = k[i] * x;
        Major[i] = sqrt(x * x + y * y);
        Minor[i] = std::abs(e[i]) / (Major[i]);

        if (Xparam[i] > 0)
            phi[i] = phi[i] + PI() / 2;
    };

    Data.resize(6);
    for (int i = 0; i < 6; i++)
        Data[i].resize(particlesPerBunch);

    double xtmp, ytmp;
    for (int i = 0; i < particlesPerBunch; i++)
    {
        if (flagsParams[0] == 0)
        {
            xtmp = (Major[0] / 3) * distribution(generator);
            ytmp = (Minor[0] / 3) * distribution(generator);
        };
        if (flagsParams[0] == 1)
        {
            while (1)
            {
                xtmp = -Major[0] + distributionUn(generator) * (2 * Major[0]);
                ytmp = -Minor[0] + distributionUn(generator) * (2 * Minor[0]);
                if ((xtmp * xtmp / (Major[0] * Major[0]) + ytmp * ytmp / (Minor[0] * Minor[0])) < 1)
                    break;
            }
        };

        if (Xparam[0] < 0)
        {
            Data[0][i] = xtmp * std::cos(phi[0]) - ytmp * std::sin(phi[0]) + DistribParams[12];
            Data[1][i] = xtmp * std::sin(phi[0]) + ytmp * std::cos(phi[0]) + DistribParams[14];
        }
        else
        {
            Data[1][i] = xtmp * std::cos(phi[0]) - ytmp * std::sin(phi[0]) + DistribParams[14];
            Data[0][i] = xtmp * std::sin(phi[0]) + ytmp * std::cos(phi[0]) + DistribParams[12];
        };

        if (flagsParams[0] == 0)
        {
            xtmp = (Major[1] / 3) * distribution(generator);
            ytmp = (Minor[1] / 3) * distribution(generator);
        };
        if (flagsParams[0] == 1)
        {
            while (1)
            {
                xtmp = -Major[1] + distributionUn(generator) * (2 * Major[0]);
                ytmp = -Minor[1] + distributionUn(generator) * (2 * Minor[0]);
                if ((xtmp * xtmp / (Major[1] * Major[1]) + ytmp * ytmp / (Minor[1] * Minor[1])) < 1)
                    break;
            }
        };

        if (Xparam[1] < 0)
        {
            Data[2][i] = xtmp * std::cos(phi[1]) - ytmp * std::sin(phi[1]) + DistribParams[13];
            Data[3][i] = xtmp * std::sin(phi[1]) + ytmp * std::cos(phi[1]) + DistribParams[15];
        }
        else
        {
            Data[3][i] = xtmp * std::cos(phi[1]) - ytmp * std::sin(phi[1]) + DistribParams[15];
            Data[2][i] = xtmp * std::sin(phi[1]) + ytmp * std::cos(phi[1]) + DistribParams[13];
        };

        if (flagsParams[0] == 0)
        {
            Data[4][i] = (DistribParams[5] + DistribParams[4]) / 2 +
                         ((DistribParams[5] - DistribParams[4]) / 6) * distribution(generator);
            Data[5][i] = P + (P * DistribParams[2] / (100 * 3)) * distribution(generator);
        };
        if (flagsParams[0] == 1)
        {
            double dp  = DistribParams[2] / 100;
            Data[4][i] = DistribParams[4] +
                         distributionUn(generator) * ((DistribParams[5] - DistribParams[4]));
            Data[5][i] = P - P * dp + distributionUn(generator) * (2 * P * dp);
        };
        while (Data[4][i] < 0)
            Data[4][i] = Data[4][i] + 2 * PI();
    }

    std::vector<int> IndexOrder(particlesPerBunch); // indexes of sorted phase
    for (int k = 0; k < particlesPerBunch; k++)
    {
        IndexOrder[k] = k;
    }
    std::sort(IndexOrder.begin(), IndexOrder.end(),
              [&](int i, int j) { return Data[4][i] > Data[4][j]; });
    std::vector<double> temp0 = Data[4];
    std::vector<double> temp1 = Data[5];
    for (int k1 = 0; k1 < particlesPerBunch; k1++) // sort energy accordingly
    {
        Data[4][k1] = temp0[IndexOrder[k1]];
        Data[5][k1] = temp1[IndexOrder[k1]];
    }
    double max = *std::max_element(Data[1].begin(), Data[1].end());
};

template <class PointType>
int EmitterDeviceBase<PointType>::GetnParticlesBunch()
{
    return NumbersParams[0];
}

template <class PointType>
const std::shared_ptr<ElectrodeCurrent<PointType>>&
EmitterDeviceBase<PointType>::GetAssignedElectrode()
{
    return Assignedelectrode;
};
template <class PointType>
void EmitterDeviceBase<PointType>::SetGetAssignedElectrode(
    const std::shared_ptr<ElectrodeCurrent<PointType>>& in)
{
    Assignedelectrode = in;
};

template <class PointType>
void EmitterDeviceBase<PointType>::SetDirectionPoints(std::vector<double> sP,
                                                      std::vector<double> eP)
{
    startPoint = sP;
    endPoint   = eP;
};

template <class PointType>
std::vector<std::vector<double>> EmitterDeviceBase<PointType>::GetDirectionPoints()
{
    std::vector<std::vector<double>> result(2);
    result[0] = startPoint;
    result[1] = endPoint;
    return result;
};

template <class PointType>
std::shared_ptr<ParticleSource2d<PointType>> EmitterDeviceBase<PointType>::GetParticleSource()
{
    return particleSource;
};
template <class PointType>
std::vector<std::shared_ptr<ParticleSource2d<PointType>>>
EmitterDeviceBase<PointType>::GetParticleSources()
{
    return std::vector<std::shared_ptr<ParticleSource2d<PointType>>>{particleSource};
};

template <class PointType>
std::vector<int> EmitterDeviceBase<PointType>::GetBoundariesList()
{
    return boundaryList;
};
template <class PointType>
EmitterDeviceBase<PointType>::EmitterDeviceBase(int DistributionStyleIn)
    : get_energy_distribution( nullptr )
{
    if (DistributionStyleIn == 0)
    {
        particleSource =
            std::shared_ptr<ParticleSource2d<PointType>>(new ParticleSource2d<PointType>());
        flagsParams.resize(1);
        NumbersParams.resize(3);
        NumbersParams[0] = 1;
        NumbersParams[1] = 1;
        NumbersParams[2] = 1;
        DistribParams.resize(2);
        DistribParams[0] = 1;
        DistribParams[1] = 1e10;
        emitterInitParam.resize(4);
        emitterInitParam[0] = 5000;
    }
    if (DistributionStyleIn == 1)
    {
        particleSource =
            std::shared_ptr<ParticleSource2d<PointType>>(new ParticleSource2d<PointType>());
        flagsParams.resize(1);
        NumbersParams.resize(3);
        NumbersParams[0] = 1;
        NumbersParams[1] = 1;
        NumbersParams[2] = 1;
        DistribParams.resize(1);
        DistribParams[0] = 1;
        emitterInitParam.resize(4);
        emitterInitParam[0] = 5000;
    }
    if (DistributionStyleIn == 5)
    {
        flagsParams.resize(1);
        NumbersParams.resize(2);
        NumbersParams[0] = 5000;
        NumbersParams[1] = 5;
        DistribParams.resize(16);
    }
    // startPoint = { 0, 0, 0 };
    // endPoint = { 0, 0, 1 };
    // startPoint.resize(3);
    // endPoint.resize(3);
};

template <class PointType>
double EmitterDeviceBase<PointType>::GetLambda()
{
    return LIGHT_VELOCITY() / DistribParams[3];
};
template <class PointType>
int EmitterDeviceBase<PointType>::GetMaxParticles()
{
    int result = 1;
    for (int i = 1; i < NumbersParams.size(); i++)
    {
        if (NumbersParams[i] > 0)
        {
            result = result * NumbersParams[i];
        };
    }
    return result;
};
template <class PointType>
int EmitterDeviceBase<PointType>::GetNumbersOfParticlesGeneration()
{
    int result = 1;
    for (int i = 1; i < NumbersParams.size(); i++)
    {
        if (NumbersParams[i] > 0)
        {
            result = result * NumbersParams[i];
        };
    }
    return result;
};
template <class PointType>
int EmitterDeviceBase<PointType>::GetEmitPeriod()
{
    return NumbersParams[0];
};

template <class PointType>
std::vector<std::vector<double>> EmitterDeviceBase<PointType>::GetAllParameters()
{
    std::vector<std::vector<double>> result(flagsParams.size() + 2);
    for (int i    = 0; i < flagsParams.size(); i++)
        result[i] = std::vector<double>{double(flagsParams[i])};

    for (int i = 0; i < NumbersParams.size(); i++)
        result[flagsParams.size()].push_back(NumbersParams[i]);

    result.back() = DistribParams;

    return result;
};
template <class PointType>
void EmitterDeviceBase<PointType>::SetAllParameters(const std::vector<std::vector<double>>& In)
{
    for (int i         = 0; i < flagsParams.size(); i++)
        flagsParams[i] = In[i][0];

    for (int i           = 0; i < NumbersParams.size(); i++)
        NumbersParams[i] = In[flagsParams.size()][i];

    DistribParams = In[flagsParams.size() + 1];
    // NumbersParams
};

template <class PointType>
std::vector<double> EmitterDeviceBase<PointType>::GetDistribParams()
{
    return DistribParams;
};

template <class PointType>
int EmitterDeviceBase<PointType>::GetEmissionType()
{
    return flagsParams[0];
};

template <class PointType>
void EmitterDeviceBase<PointType>::SetFlowCurrent(double res)
{
    particleSource->SetFlowCurrent(res);
};
template <class PointType>
double EmitterDeviceBase<PointType>::GetSourceSize()
{
    if (particleSource->sourceSurface.size() == 0)
        return 0;
    return particleSource->sourceSurface.back().curveLength;
};
template <class PointType>
int EmitterDeviceBase<PointType>::GetnParticlesEmitter()
{
    return NumbersParams[2];
};
template <class PointType>
double EmitterDeviceBase<PointType>::getErAverage()
{
    return particleSource->getErAverage();
};
template <class PointType>
std::vector<std::vector<double>> EmitterDeviceBase<PointType>::GetEmitterField()
{
    return particleSource->GetEmitterField();
};
template <class PointType>
std::vector<double> EmitterDeviceBase<PointType>::GetEmitterInitParameters()
{
    return emitterInitParam;
};

template void
EmitterDeviceBase<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void EmitterDeviceBase<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void EmitterDeviceBase<double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void EmitterDeviceBase<float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void EmitterDeviceBase<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& NumbersParams;
    ar& DistribParams;
    ar& boundaryList;
    ar& particleSource;
    ar& flagsParams;
    ar& Assignedelectrode;
    ar& emitterInitParam;
};
template <class PointType>
template <class Archive>
void EmitterDeviceBase<PointType>::load(Archive& ar, const unsigned int)
{
    ar& NumbersParams;
    ar& DistribParams;
    ar& boundaryList;
    ar& particleSource;
    ar& flagsParams;
    ar& Assignedelectrode;
    ar& emitterInitParam;
};