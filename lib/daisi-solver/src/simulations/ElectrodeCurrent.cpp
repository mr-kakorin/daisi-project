#include "ElectrodeCurrent.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Geom.h"
#include "geomTools.h"
#include <Constants.h>

template class ElectrodeCurrent<float>;
template class ElectrodeCurrent<double>;

template <class PointType>
std::vector<double> ElectrodeCurrent<PointType>::GetParameteres()
{
    return parameters;
};

template <class PointType>
void ElectrodeCurrent<PointType>::SetParameteres(const std::vector<double>& input)
{
    parameters = input;
};

template <class PointType>
ElectrodeCurrent<PointType>::ElectrodeCurrent()
{
    parameters.resize(2);
    parameters[0] = 200;
    parameters[1] = 20;
};

template <class PointType>
double ElectrodeCurrent<PointType>::GetCurrent()
{
    return totalCurrent;
};
template <class PointType>
void ElectrodeCurrent<PointType>::SetCurrent(PointType cur)
{
    totalCurrent = cur;
};
template <class PointType>
void ElectrodeCurrent<PointType>::InitParallel(int numThreads)
{
    mapSize = 10;
    /*std::vector<PointType> tmp(averageIrradiatedPowerDensity.size());
    for (int i = 0; i < mapSize; i++)
    {
            CurrentAtPointsTimes.push_back(tmp);
            EnergyAtPointsTimes.push_back(tmp);
    }*/
    EnergyAtPointsThreaded.resize(numThreads);
    CurrentAtPointsThreaded.resize(numThreads);

    CurrentAtPointsAverage.resize(averageIrradiatedPowerDensity.size());
    EnergyAtPointsAverage.resize(averageIrradiatedPowerDensity.size());

    for (int i = 0; i < numThreads; i++)
    {
        EnergyAtPointsThreaded[i].resize(averageIrradiatedPowerDensity.size());
        CurrentAtPointsThreaded[i].resize(averageIrradiatedPowerDensity.size());
        for (int j = 0; j < EnergyAtPointsThreaded[i].size(); j++)
        {
            EnergyAtPointsThreaded[i][j]  = 0;
            CurrentAtPointsThreaded[i][j] = 0;
        };
    }

    for (int j = 0; j < averageIrradiatedPowerDensity.size(); j++)
    {
        averageIrradiatedPowerDensity[j]   = 0;
        averageIrradiatedCurrentDensity[j] = 0;
        averageCollectedCurrentDensity[j]  = 0;
        CurrentAtPointsAverage[j]          = 0;
        EnergyAtPointsAverage[j]           = 0;
    }
    CurrentAtPointsLongTimes.clear();
    EnergyAtPointsLongTimes.clear();
    CurrentAtPointsTimes.clear();
    EnergyAtPointsTimes.clear();
};
template <class PointType>
void ElectrodeCurrent<PointType>::SetBoundaryList(
    std::vector<int> list, std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
    std::string& errorMsg){

};

/*template std::vector<std::vector<float>>  ElectrodeCurrent<double>::GetElectrodeValue<float>(int
flag); template std::vector<std::vector<double>>
ElectrodeCurrent<double>::GetElectrodeValue<double>(int flag); template
std::vector<std::vector<float>>  ElectrodeCurrent<float>::GetElectrodeValue<float>(int flag);
template std::vector<std::vector<double>> ElectrodeCurrent<float>::GetElectrodeValue<double>(int
flag);

template <class PointType>
template <class PointType1>*/

template <class PointType>
std::vector<std::vector<double>> ElectrodeCurrent<PointType>::GetElectrodeValueD(int flag)
{
    std::vector<std::vector<double>> result;
    result.resize(2);

    if (ElectrodeEdges.size() == 0)
        return result;

    result[0].resize(ElectrodeEdges.size() - 1);
    result[1].resize(ElectrodeEdges.size() - 1);

    if (ElectrodeEdges.size() > 0)
    {
        float l = 0;
        // result[0][0] = 0;
        for (int i = 1; i < ElectrodeEdges.size(); i++)
        {
            result[0][i - 1] = l;
            l                = l + ElectrodeEdges[i - 1].length();
            if (flag == 0)
                result[1][i - 1] = averageIrradiatedPowerDensity[i];
            if (flag == 1)
                result[1][i - 1] = averageIrradiatedCurrentDensity[i];
            if (flag == 2)
                result[1][i - 1] = averageCollectedCurrentDensity[i];
        };
    }
    return result;
};
template <class PointType>
std::vector<std::vector<float>> ElectrodeCurrent<PointType>::GetElectrodeValueF(int flag)
{
    std::vector<std::vector<float>> result;
    result.resize(2);

    if (ElectrodeEdges.size() == 0)
        return result;

    result[0].resize(ElectrodeEdges.size() - 1);
    result[1].resize(ElectrodeEdges.size() - 1);

    if (ElectrodeEdges.size() > 0)
    {
        float l = 0;
        // result[0][0] = 0;
        for (int i = 1; i < ElectrodeEdges.size(); i++)
        {
            result[0][i - 1] = l;
            l                = l + ElectrodeEdges[i - 1].length();
            if (flag == 0)
                result[1][i - 1] = averageIrradiatedPowerDensity[i];
            if (flag == 1)
                result[1][i - 1] = averageIrradiatedCurrentDensity[i];
            if (flag == 2)
                result[1][i - 1] = averageCollectedCurrentDensity[i];
        };
    }
    return result;
};

template <class PointType>
void ElectrodeCurrent<PointType>::AddCharge(PointType charge, const DGeo::Point<PointType>& point,
                                            PointType dt, PointType energy,
                                            PointType chargeParticle, int thread)
{
    if (std::isnan(charge) || std::isinf(charge) || std::isnan(dt) || std::isinf(dt))
        return;

    PointType tol     = 1e-9;
    PointType minDist = 1e9;
    int       i       = 0;
    for (int j = 0; j < ElectrodeEdges.size() - 1; j++)
    {

        PointType d1 = point.Dist2Point(ElectrodeEdges[j].point1);
        PointType d2 = point.Dist2Point(ElectrodeEdges[j].point2);

        if (d1 + d2 < minDist)
        {
            minDist = d1;
            i       = j;
        }
    }

    PointType dist = point.Dist2Point(ElectrodeEdges[i].point1);

    PointType d1;
    PointType d2;

    PointType l  = ElectrodeEdges[i].length();
    PointType w1 = (l - dist) / l;
    if (std::abs(w1) > 1 || std::isnan(w1) || std::isinf(w1))
        w1 = 1;

    PointType w2 = 1 - w1;

    CurrentAtPointsThreaded[thread][i] = CurrentAtPointsThreaded[thread][i] + w1 * charge / dt;
    EnergyAtPointsThreaded[thread][i] =
        EnergyAtPointsThreaded[thread][i] + std::abs(w1 * (charge * energy / (dt * chargeParticle)));

    CurrentAtPointsThreaded[thread][i + 1] =
        CurrentAtPointsThreaded[thread][i + 1] + w2 * charge / dt;
    EnergyAtPointsThreaded[thread][i + 1] =
        EnergyAtPointsThreaded[thread][i + 1] + std::abs(w2 * (charge * energy / (dt * chargeParticle)));

    /*if (i == 0 || d1 < d2)
    {
            PointType l = ElectrodeEdges[i].length();
            PointType w1 = (l - dist) / l;
            if (std::abs(w1) > 1 || std::isnan(w1) || std::isinf(w1))
                    w1 = 1;

            PointType w2 = 1 - w1;






            CurrentAtPointsThreaded[thread][i + 1] = CurrentAtPointsThreaded[thread][i + 1] +
    w2*charge / dt; EnergyAtPointsThreaded[thread][i + 1] = EnergyAtPointsThreaded[thread][i + 1] +
    std::abs(w2*(charge * energy / (dt *chargeParticle)));


            CurrentAtPointsThreaded[thread][i] = CurrentAtPointsThreaded[thread][i] + w1*charge /
    dt; EnergyAtPointsThreaded[thread][i] = EnergyAtPointsThreaded[thread][i] + std::abs(w1*(charge *
    energy / (dt *chargeParticle))) ;

    }
    else
    {
            PointType l = ElectrodeEdges[i-1].length();

            PointType w1 = (l - dist) / l;
            if (std::abs(w1) > 1 || std::isnan(w1) || std::isinf(w1))
                    w1 = 1;

            PointType w2 = 1 - w1;

            CurrentAtPointsThreaded[thread][i] = CurrentAtPointsThreaded[thread][i] + w1*(charge) /
    dt; EnergyAtPointsThreaded[thread][i] = EnergyAtPointsThreaded[thread][i] + std::abs(w1*(charge *
    energy / (dt *chargeParticle)));


            CurrentAtPointsThreaded[thread][i - 1] = CurrentAtPointsThreaded[thread][i - 1] +
    w2*charge / dt; EnergyAtPointsThreaded[thread][i - 1] = EnergyAtPointsThreaded[thread][i - 1] +
    std::abs(w2*(charge  * energy / (dt * chargeParticle)));

    };*/
}
template <class PointType>
void ElectrodeCurrent<PointType>::ResetCurrent()
{
    /*for (int i = 0; i < CurrentAtPoints.size(); i++)
    {
            for (int j = 0; j < CurrentAtPoints[i].size(); j++)
            {
                    CurrentAtPoints[i][j] = 0;
                    CurrentDensityAtPoints[i][j] = 0;
            }
    }*/
}

template <class PointType>
void ElectrodeCurrent<PointType>::ResetPower()
{
    CurrentAtPointsLongTimes.clear();
    EnergyAtPointsLongTimes.clear();
    CurrentAtPointsTimes.clear();
    EnergyAtPointsTimes.clear();
    for (int i = 0; i < CurrentAtPointsThreaded.size(); i++)
    {
        for (int j = 0; j < EnergyAtPointsThreaded[i].size(); j++)
        {
            EnergyAtPointsThreaded[i][j]  = 0;
            CurrentAtPointsThreaded[i][j] = 0;
        };
    }

    for (int j = 0; j < averageIrradiatedPowerDensity.size(); j++)
    {
        averageIrradiatedPowerDensity[j]   = 0;
        averageIrradiatedCurrentDensity[j] = 0;
        averageCollectedCurrentDensity[j]  = 0;
        CurrentAtPointsAverage[j]          = 0;
        EnergyAtPointsAverage[j]           = 0;
    }
}

template <class PointType>
void ElectrodeCurrent<PointType>::PowerAndCurrentsCalculate(int k, double t, double dt, int type)
{
    if (k == -1)
        return;

    int n = parameters[1] / dt;
    for (int i = 1; i < CurrentAtPointsThreaded.size(); i++)
    {
        for (int j = 0; j < CurrentAtPointsThreaded[i].size(); j++)
        {
            CurrentAtPointsThreaded[0][j] =
                CurrentAtPointsThreaded[0][j] + CurrentAtPointsThreaded[i][j];
            CurrentAtPointsThreaded[i][j] = 0;

            EnergyAtPointsThreaded[0][j] =
                EnergyAtPointsThreaded[0][j] + EnergyAtPointsThreaded[i][j];
            EnergyAtPointsThreaded[i][j] = 0;
        }
    }

    if (mapSize <= CurrentAtPointsTimes.size())
    {
        CurrentAtPointsTimes.pop_front();
        EnergyAtPointsTimes.pop_front();
    }
    CurrentAtPointsTimes.push_back(CurrentAtPointsThreaded[0]);
    EnergyAtPointsTimes.push_back(EnergyAtPointsThreaded[0]);

    double steps = std::min(k + 1, n);

    int sizeLong = n / mapSize;

    for (int j = 0; j < averageIrradiatedCurrentDensity.size(); j++)
    {
        CurrentAtPointsAverage[j] = 0;
        EnergyAtPointsAverage[j]  = 0;
        for (int m = 0; m < CurrentAtPointsTimes.size(); m++)
        {
            CurrentAtPointsAverage[j] = CurrentAtPointsAverage[j] + CurrentAtPointsTimes[m][j];
            EnergyAtPointsAverage[j]  = EnergyAtPointsAverage[j] + EnergyAtPointsTimes[m][j];
        }
        CurrentAtPointsAverage[j]     = CurrentAtPointsAverage[j] / CurrentAtPointsTimes.size();
        CurrentAtPointsThreaded[0][j] = 0;
        EnergyAtPointsAverage[j]      = EnergyAtPointsAverage[j] / CurrentAtPointsTimes.size();
        EnergyAtPointsThreaded[0][j]  = 0;
    };

    CurrentAtPointsAverageCollected = CurrentAtPointsAverage;

    for (int j = 1; j < CurrentAtPointsAverage.size(); j++)
        CurrentAtPointsAverageCollected[j] =
            CurrentAtPointsAverageCollected[j] + CurrentAtPointsAverageCollected[j - 1];

    averageCollectedCurrentDensitySim = CurrentAtPointsAverageCollected;

    for (int j = 0; j < CurrentAtPointsThreaded[0].size(); j++)
    {
        PointType koef;
        if (type == 2)
            koef = 2 * PI() * ElectrodeEdges[j].point1.x;
        else
            koef                             = 1;
        averageCollectedCurrentDensitySim[j] = averageCollectedCurrentDensitySim[j] / koef;
    };

    SetCurrent(std::abs(CurrentAtPointsAverageCollected.back()));

    if (0 == k % mapSize)
    {

        if (sizeLong <= CurrentAtPointsLongTimes.size() && !CurrentAtPointsLongTimes.empty())
        {
            CurrentAtPointsLongTimes.pop_front();
            EnergyAtPointsLongTimes.pop_front();
        }

        CurrentAtPointsLongTimes.push_back(CurrentAtPointsAverage);
        EnergyAtPointsLongTimes.push_back(EnergyAtPointsAverage);

        for (int j = 0; j < averageIrradiatedCurrentDensity.size(); j++)
        {
            averageIrradiatedCurrentDensity[j] = 0;
            averageIrradiatedPowerDensity[j]   = 0;
            for (int m = 0; m < EnergyAtPointsLongTimes.size(); m++)
            {
                averageIrradiatedCurrentDensity[j] =
                    averageIrradiatedCurrentDensity[j] + CurrentAtPointsLongTimes[m][j];
                averageIrradiatedPowerDensity[j] =
                    averageIrradiatedPowerDensity[j] + EnergyAtPointsLongTimes[m][j];
            }
            averageIrradiatedCurrentDensity[j] =
                averageIrradiatedCurrentDensity[j] / CurrentAtPointsLongTimes.size();
            averageIrradiatedPowerDensity[j] =
                averageIrradiatedPowerDensity[j] / EnergyAtPointsLongTimes.size();
        }
        averageCollectedCurrentDensity = averageIrradiatedCurrentDensity;

        for (int j = 1; j < CurrentAtPointsAverage.size(); j++)
            averageCollectedCurrentDensity[j] =
                averageCollectedCurrentDensity[j] + averageCollectedCurrentDensity[j - 1];

        for (int j = 0; j < averageIrradiatedCurrentDensity.size(); j++)
        {
            PointType koef;
            if (type == 2)
                koef = 2 * PI() * ElectrodeEdges[j].point1.x;
            else
                koef = 1;
            averageIrradiatedCurrentDensity[j] =
                std::abs(averageIrradiatedCurrentDensity[j] / (ElectrodeEdges[j].length() * koef));
            averageIrradiatedPowerDensity[j] =
                (averageIrradiatedPowerDensity[j] / (ElectrodeEdges[j].length() * koef));
            averageCollectedCurrentDensity[j] = averageCollectedCurrentDensity[j] / koef;
        };

        averageIrradiatedPowerDensity[0]   = averageIrradiatedPowerDensity[1];
        averageIrradiatedCurrentDensity[0] = averageIrradiatedCurrentDensity[1];
    };
}

template <class PointType>
void ElectrodeCurrent<PointType>::SetBoundaryList(
    std::vector<int> list, std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
    std::string& errorMsg)
{
    boundaryNumbers = list;

    ElectrodeEdges.clear();
    CurrentAtPointsThreaded.resize(1);
    EnergyAtPointsThreaded.resize(1);

    if (!list.size())
        return;

    mergeSortResize(parameters[0], list, boundaries, ElectrodeEdges, errorMsg);

    if (errorMsg.size())
        return;

    averageIrradiatedPowerDensity.resize(ElectrodeEdges.size());
    averageIrradiatedCurrentDensity.resize(ElectrodeEdges.size());
    averageCollectedCurrentDensity.resize(ElectrodeEdges.size());

    /*for (int i = 0; i < ElectrodeEdges.size(); i++)
    {
            ElectrodeEdges[i].point1 = ElectrodeEdges[i].Middle();
    };*/

    /*int normSize = 5000;

    std::vector<DGeo::Edge<PointType>> EdgesData = boundary.EdgesData;

    if (EdgesData.size()>1)
            std::sort(EdgesData.begin(), EdgesData.end(), EdgeCmp);

    //if (size >= normSize)
    double length = 0;
    for (int i = 0; i < size; i++)
    {
            length = length + EdgesData[i].length();
            //	currentDensityDistribution.push_back();
    };
    double avLength = length / normSize;
    double curveLength = 0;
    for (int i = 0; i < size; i++)
    {
            if (EdgesData[i].length() < avLength)
            {
                    sourceSurface.push_back(*(new ParticleSourceEdge<PointType>(EdgesData[i],
    curveLength))); curveLength = curveLength + EdgesData[i].length();
            }
            else
            {
                    int resize = ceil(EdgesData[i].length() / avLength);
                    std::vector<DGeo::Edge <PointType>> resArr = EdgesData[i].resize(resize);
                    for (int j = 0; j < resize; j++)
                    {
                            sourceSurface.push_back(*(new ParticleSourceEdge<PointType>(resArr[j],
    curveLength))); curveLength = curveLength + resArr[j].length();
                    };
            };
    };*/
};

template void
ElectrodeCurrent<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);
template void ElectrodeCurrent<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void
ElectrodeCurrent<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void ElectrodeCurrent<float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void ElectrodeCurrent<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& boundaryNumbers;
    ar& ElectrodeEdges;
    ar& parameters;
    ar& averageIrradiatedPowerDensity;
    ar& averageIrradiatedCurrentDensity;
    ar& averageCollectedCurrentDensity;
}
template <class PointType>
template <class Archive>
void ElectrodeCurrent<PointType>::load(Archive& ar, const unsigned int)
{

    ar& boundaryNumbers;
    ar& ElectrodeEdges;
    ar& parameters;
    ar& averageIrradiatedPowerDensity;
    ar& averageIrradiatedCurrentDensity;
    ar& averageCollectedCurrentDensity;
    /*
    CurrentAtPoints.resize(1);
    parameters.resize(2);
    CurrentDensityAtPoints.resize(1);
    ar & boundaryNumbers;
    ar & ElectrodeEdges;
    ar & CurrentAtPoints[0];
    ar & PowerAtPoints;
    ar & CurrentDensityAtPoints[0];
    PowerAtPointsTmp.resize(1);
    EnergyAtPoints.resize(1);
    EnergyAtPoints[0].resize(PowerAtPoints.size());
    PowerAtPointsAverage.resize(PowerAtPoints.size());

    averageIrradiatedPowerDensit.resize(PowerAtPoints.size());
    averageIrradiatedCurrentDensity.resize(PowerAtPoints.size());
    averageCollectedCurrentDensity.resize(PowerAtPoints.size());

    for (int j = 0; j < EnergyAtPoints[0].size(); j++)
            EnergyAtPoints[0][j] = 0;*/
}